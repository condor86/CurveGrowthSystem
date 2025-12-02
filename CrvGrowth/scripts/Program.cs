// File: CrvGrowth/Program.cs
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Numerics;

using NSGAII;
using CrvGrowth;
using CrvGrowth.Scripts;

namespace CrvGrowth
{
    class Program
    {
        private enum OptimizationMode
        {
            CrvGrowth,
            Blinds
        }

        // 在这里切换这次要跑的优化类型：
        //private const OptimizationMode Mode = OptimizationMode.CrvGrowth;
        private const OptimizationMode Mode = OptimizationMode.Blinds;

        // —— 站点参数（与 LightingSimulator 默认一致：南京）——
        private const double SiteLatitudeDeg   = 32.0603;
        private const double SiteLongitudeDeg  = 118.7969;
        private const double SiteTimezoneHours = 8.0;

        // —— 模型坐标系：Up=+Z, North=+Y —— 
        private static readonly Vector3 Up    = new(0, 0, 1);
        private static readonly Vector3 North = new(0, 1, 0);

        static void Main(string[] args)
        {
            var totalWatch = Stopwatch.StartNew();

            // 目录结构
            string rootDir   = AppDomain.CurrentDomain.BaseDirectory;
            string parentDir = Path.GetFullPath(Path.Combine(rootDir, "..", "..", ".."));
            string dataDir   = Path.Combine(parentDir, "data");
            string resultDir = Path.Combine(parentDir, "results");
            Directory.CreateDirectory(resultDir);

            switch (Mode)
            {
                case OptimizationMode.CrvGrowth:
                    RunCrvGrowthOptimization(dataDir, resultDir);
                    break;
                case OptimizationMode.Blinds:
                    RunBlindsOptimization(resultDir);
                    break;
                default:
                    throw new InvalidOperationException("未知的优化模式。");
            }

            totalWatch.Stop();
            Console.WriteLine($"All done. Total time: {totalWatch.Elapsed}");
        }

        // =====================================================================
        // 一、原有 CrvGrowth 优化流程（保持原逻辑）
        // =====================================================================

        private static void RunCrvGrowthOptimization(string dataDir, string resultDir)
        {
            // 输入数据路径
            string startingCsv = Path.Combine(dataDir, "iStartingPositions.csv");
            string repellerCsv = Path.Combine(dataDir, "iRepellers.csv");

            // 预加载输入
            var startingPoints = IOHelper.LoadPointsFromFile(startingCsv);
            var repellerPoints = IOHelper.LoadPointsFromFile(repellerCsv);

            // 预计算太阳向量（一次）：夏/冬各一组
            var summerToSuns = SunVectors.Build(
                NSGAWiring.SummerDate, NSGAWiring.StartTime, NSGAWiring.EndTime, NSGAWiring.Interval,
                SiteLatitudeDeg, SiteLongitudeDeg, SiteTimezoneHours,
                Up, North);

            var winterToSuns = SunVectors.Build(
                NSGAWiring.WinterDate, NSGAWiring.StartTime, NSGAWiring.EndTime, NSGAWiring.Interval,
                SiteLatitudeDeg, SiteLongitudeDeg, SiteTimezoneHours,
                Up, North);

            // NSGA-II 基因边界（4 + 400 + 400 = 804）
            const int repellerCount = 4;
            const int offsetCount   = 400;
            const int angleCount    = 400;
            const int geneLen       = repellerCount + offsetCount + angleCount;

            var lo = new double[geneLen];
            var hi = new double[geneLen];

            // 1) 4 个 repeller 因子
            for (int i = 0; i < repellerCount; i++)
            {
                lo[i] = 0.01;
                hi[i] = 5.0;
            }

            // 2) 400 个逐点 -Y 偏移
            for (int i = repellerCount; i < repellerCount + offsetCount; i++)
            {
                lo[i] = 50.0;
                hi[i] = 100.0;
            }

            // 3) 400 个逐点旋转角（度）
            for (int i = repellerCount + offsetCount; i < geneLen; i++)
            {
                lo[i] = -25.0;
                hi[i] = 25.0;
            }

            // NSGA-II 日志目录
            string nsgaLogDir = Path.Combine(resultDir, "nsga_logs_crvgrowth");
            Directory.CreateDirectory(nsgaLogDir);

            // NSGA-II 配置（向量路径）
            var cfg = new NSGAConfig
            {
                PopulationSize = 50,   // 正式跑可以改 50
                Generations    = 100,   // 正式跑可以改 100
                CrossoverRate  = 0.9,
                MutationRate   = 1.0 / geneLen,
                GeneLength     = geneLen,
                LowerBounds    = lo,
                UpperBounds    = hi,
                RandomSeed     = 1,
                DegreeOfParallelism = Environment.ProcessorCount,
                LogDir         = nsgaLogDir,
                SbxEta         = 20.0,
                PolyMutationEta= 20.0,
                Evaluate       = NSGAWiring.MakeEvaluator(
                    startingPoints, repellerPoints,
                    summerToSuns, winterToSuns)
            };

            Console.WriteLine("NSGA-II optimization (CrvGrowth) started...");
            var runWatch = Stopwatch.StartNew();

            var solver   = new NSGAII.NSGAII(cfg);
            var finalPop = solver.Run();

            runWatch.Stop();
            Console.WriteLine($"NSGA-II (CrvGrowth) finished in {runWatch.Elapsed}.");

            // 最终一代 Pareto 前沿 & 代表解
            var pareto = finalPop.Where(ind => ind.Rank == 0).ToList();
            Console.WriteLine($"[CrvGrowth] Final Pareto size: {pareto.Count}");

            var rep = pareto.OrderBy(ind => ind.Objectives.Sum()).First();
            Console.WriteLine("[CrvGrowth] Exporting representative solution geometry & lighting...");

            // 导出 CrvGrowth 解
            string outCrvCsv              = Path.Combine(resultDir, "resultsCrv.csv");
            string outVerticalCsv         = Path.Combine(resultDir, "resultsVertical.csv");
            string outFilletVerticalCsv   = Path.Combine(resultDir, "resultsFilletVertical.csv");
            string outLightingSummer      = Path.Combine(resultDir, "resultsLighting_summer.csv");
            string outLightingWinter      = Path.Combine(resultDir, "resultsLighting_winter.csv");
            string outNurbsCsv            = Path.Combine(resultDir, "resultsNurbs.csv");
            string outFilletCsv           = Path.Combine(resultDir, "resultsFillet.csv");

            SaveSolutionGeometryAndLighting(
                genes: rep.Genes,
                startingPoints: startingPoints,
                repellerPoints: repellerPoints,
                outVerticalCsv: outVerticalCsv,
                outFilletVerticalCsv: outFilletVerticalCsv,
                outCrvCsv: outCrvCsv,
                outLightingSummerCsv: outLightingSummer,
                outLightingWinterCsv: outLightingWinter,
                outNurbsCsv: outNurbsCsv,
                outFilletCsv: outFilletCsv,
                summerToSuns: summerToSuns,
                winterToSuns: winterToSuns
            );
        }

        // CrvGrowth：导出几何与光照
        private static void SaveSolutionGeometryAndLighting(
            double[] genes,
            List<Vector3> startingPoints,
            List<Vector3> repellerPoints,
            string outVerticalCsv,
            string outFilletVerticalCsv,
            string outCrvCsv,
            string outLightingSummerCsv,
            string outLightingWinterCsv,
            string outNurbsCsv,
            string outFilletCsv,
            Vector3[] summerToSuns,
            Vector3[] winterToSuns)
        {
            const int repellerCount = 4;
            const int offsetCount   = 400;
            const int angleCount    = 400;
            const float filletRadiusDefault = 10f;

            var repellerFactors = genes.Take(repellerCount).ToList();
            var offsets         = genes.Skip(repellerCount).Take(offsetCount).ToArray();
            var anglesDeg       = genes.Skip(repellerCount + offsetCount).Take(angleCount).ToArray();

            // 生长（平面）
            var system = new GrowthSystem();
            var flatCurve = system.Run(
                starting:        startingPoints,
                repellers:       repellerPoints,
                repellerFactors: repellerFactors,
                maxPointCount:   NSGAWiring.MaxPointCount,
                maxIterCount:    NSGAWiring.MaxIterCount,
                baseDist:        NSGAWiring.BaseDist
            );

            // 转垂直 (x, y, 0) → (x, 0, z=y)
            var verticalCrv = flatCurve.Select(p => new Vector3(p.X, 0f, p.Y)).ToList();
            IOHelper.SavePointsToFile(outVerticalCsv, verticalCrv);

            // 垂直曲线圆角化
            var filletVertical = FilletUtil.FilletPolylineWithFixedArcPoints(
                pts: verticalCrv,
                radius: filletRadiusDefault,
                arcPointCount: 9,
                angleEpsDeg: 1.0f,
                isClosed: true,
                clampRadius: true);
            IOHelper.SavePointsToFile(outFilletVerticalCsv, filletVertical);

            // 逐点 -Y 偏移
            int N = Math.Min(verticalCrv.Count, offsets.Length);
            var extrudedCrv = new List<Vector3>(verticalCrv);
            for (int i = 0; i < N; i++)
            {
                var p = verticalCrv[i];
                extrudedCrv[i] = new Vector3(p.X, p.Y - (float)offsets[i], p.Z);
            }

            // 局部平面旋转
            NSGAWiring.ApplyLocalPlaneRotation(verticalCrv, extrudedCrv, anglesDeg);

            // 挤出后曲线
            IOHelper.SavePointsToFile(outCrvCsv, extrudedCrv);

            // 夏/冬光照（向量直跑）
            SimAndSaveVectors(verticalCrv, extrudedCrv, summerToSuns, outLightingSummerCsv);
            SimAndSaveVectors(verticalCrv, extrudedCrv, winterToSuns, outLightingWinterCsv);

            Console.WriteLine($"Saved: {outVerticalCsv}");
            Console.WriteLine($"Saved: {outFilletVerticalCsv}");
            Console.WriteLine($"Saved: {outCrvCsv}");
            Console.WriteLine($"Saved: {outLightingSummerCsv}");
            Console.WriteLine($"Saved: {outLightingWinterCsv}");

            var filletSampled = FilletUtil.FilletPolylineWithFixedArcPoints(
                pts: extrudedCrv,
                radius: filletRadiusDefault,
                arcPointCount: 9,
                angleEpsDeg: 1.0f,
                isClosed: true,
                clampRadius: true);

            IOHelper.SavePointsToFile(outFilletCsv, filletSampled);
            Console.WriteLine($"Saved (fillet points): {outFilletCsv}");
        }

        private static void SimAndSaveVectors(
            List<Vector3> verticalCrv,
            List<Vector3> extrudedCrv,
            Vector3[] toSuns,
            string outCsv)
        {
            var sim = new LightingSimulator(
                verticalCurve: verticalCrv,
                extrudedCurve: extrudedCrv,
                date:          NSGAWiring.SummerDate,
                startTime:     NSGAWiring.StartTime,
                endTime:       NSGAWiring.EndTime,
                interval:      NSGAWiring.Interval,
                roomWidth:     NSGAWiring.RoomWidth,
                roomDepth:     NSGAWiring.RoomDepth,
                gridSize:      NSGAWiring.GridSize,
                isClosed:      true
            );

            sim.RunWithSunVectors(toSuns);
            sim.SaveLightHourGrid(outCsv);
        }

        // =====================================================================
        // 二、新百叶窗优化流程
        // =====================================================================

        private static void RunBlindsOptimization(string resultDir)
        {
            // 预计算太阳向量（夏/冬）
            var summerToSuns = SunVectors.Build(
                NSGAWiringBlinds.SummerDate, NSGAWiringBlinds.StartTime, NSGAWiringBlinds.EndTime, NSGAWiringBlinds.Interval,
                SiteLatitudeDeg, SiteLongitudeDeg, SiteTimezoneHours,
                Up, North);

            var winterToSuns = SunVectors.Build(
                NSGAWiringBlinds.WinterDate, NSGAWiringBlinds.StartTime, NSGAWiringBlinds.EndTime, NSGAWiringBlinds.Interval,
                SiteLatitudeDeg, SiteLongitudeDeg, SiteTimezoneHours,
                Up, North);

            // 基因参数（1 + 100 + 100 = 201）
            int maxBlade = NSGAWiringBlinds.MaxBladeCount;
            int geneLen  = 1 + maxBlade + maxBlade; // [0]=bladeCount, [1..100]=lengths, [101..200]=angles

            var lo = new double[geneLen];
            var hi = new double[geneLen];

            // 0) 叶片个数
            lo[0] = 2.0;
            hi[0] = maxBlade;

            // 1..maxBlade：拉伸长度（原始值，内部会再按 bladeCount clamp）
            for (int i = 1; i <= maxBlade; i++)
            {
                lo[i] = 0.0;
                hi[i] = NSGAWiringBlinds.UnitSize;
            }

            // 1+maxBlade..geneLen-1：角度
            for (int i = 1 + maxBlade; i < geneLen; i++)
            {
                lo[i] = -85.0;
                hi[i] = 85.0;
            }

            string nsgaLogDir = Path.Combine(resultDir, "nsga_logs_blinds");
            Directory.CreateDirectory(nsgaLogDir);

            var cfg = new NSGAConfig
            {
                PopulationSize = 50,   // 正式跑可改 50
                Generations    = 100,   // 正式跑可改 100
                CrossoverRate  = 0.9,
                MutationRate   = 1.0 / geneLen,
                GeneLength     = geneLen,
                LowerBounds    = lo,
                UpperBounds    = hi,
                RandomSeed     = 1,
                DegreeOfParallelism = Environment.ProcessorCount,
                LogDir         = nsgaLogDir,
                SbxEta         = 20.0,
                PolyMutationEta= 20.0,
                Evaluate       = NSGAWiringBlinds.MakeEvaluatorUsingVectors(
                    summerToSuns, winterToSuns)
            };

            Console.WriteLine("NSGA-II optimization (Blinds) started...");
            var runWatch = Stopwatch.StartNew();

            var solver   = new NSGAII.NSGAII(cfg);
            var finalPop = solver.Run();

            runWatch.Stop();
            Console.WriteLine($"NSGA-II (Blinds) finished in {runWatch.Elapsed}.");

            var pareto = finalPop.Where(ind => ind.Rank == 0).ToList();
            Console.WriteLine($"[Blinds] Final Pareto size: {pareto.Count}");

            var rep = pareto.OrderBy(ind => ind.Objectives.Sum()).First();
            Console.WriteLine("[Blinds] Exporting representative blinds geometry & lighting...");

            string outVerticalCsv    = Path.Combine(resultDir, "blinds_vertical.csv");
            string outExtrudedCsv    = Path.Combine(resultDir, "blinds_extruded.csv");
            string outLightingSummer = Path.Combine(resultDir, "blinds_lighting_summer.csv");
            string outLightingWinter = Path.Combine(resultDir, "blinds_lighting_winter.csv");

            SaveBlindsSolution(
                genes: rep.Genes,
                outVerticalCsv: outVerticalCsv,
                outExtrudedCsv: outExtrudedCsv,
                outLightingSummerCsv: outLightingSummer,
                outLightingWinterCsv: outLightingWinter,
                summerToSuns: summerToSuns,
                winterToSuns: winterToSuns
            );
        }

        private static void SaveBlindsSolution(
            double[] genes,
            string outVerticalCsv,
            string outExtrudedCsv,
            string outLightingSummerCsv,
            string outLightingWinterCsv,
            Vector3[] summerToSuns,
            Vector3[] winterToSuns)
        {
            int maxBlade = NSGAWiringBlinds.MaxBladeCount;

            // 解码 bladeCount（与 NSGAWiringBlinds 内部保持一致）
            int rawCount   = (int)Math.Round(genes[0]);
            int bladeCount = rawCount;
            if (bladeCount < 2) bladeCount = 2;
            if (bladeCount > maxBlade) bladeCount = maxBlade;

            BlindsGenerator.GenerateFromGenes(
                bladeCount:       bladeCount,
                unitSize:         NSGAWiringBlinds.UnitSize,
                genes:            genes,
                lengthStartIndex: 1,
                angleStartIndex:  1 + maxBlade,
                out var basePoints,
                out var extrudedEdge
            );

            IOHelper.SavePointsToFile(outVerticalCsv, basePoints);
            IOHelper.SavePointsToFile(outExtrudedCsv, extrudedEdge);

            SimAndSaveBlindsWithVectors(basePoints, extrudedEdge, summerToSuns, outLightingSummerCsv);
            SimAndSaveBlindsWithVectors(basePoints, extrudedEdge, winterToSuns, outLightingWinterCsv);

            Console.WriteLine($"Saved [Blinds]: {outVerticalCsv}");
            Console.WriteLine($"Saved [Blinds]: {outExtrudedCsv}");
            Console.WriteLine($"Saved [Blinds]: {outLightingSummerCsv}");
            Console.WriteLine($"Saved [Blinds]: {outLightingWinterCsv}");
        }

        private static void SimAndSaveBlindsWithVectors(
            List<Vector3> verticalCrv,
            List<Vector3> extrudedCrv,
            Vector3[] toSuns,
            string outCsv)
        {
            var sim = new LightingSimulator(
                verticalCurve: verticalCrv,
                extrudedCurve: extrudedCrv,
                date:          NSGAWiringBlinds.SummerDate,  // 占位
                startTime:     NSGAWiringBlinds.StartTime,
                endTime:       NSGAWiringBlinds.EndTime,
                interval:      NSGAWiringBlinds.Interval,
                roomWidth:     NSGAWiringBlinds.RoomWidth,
                roomDepth:     NSGAWiringBlinds.RoomDepth,
                gridSize:      NSGAWiringBlinds.GridSize,
                isClosed:      false,
                enablePeriodicTiling: true,
                mirrorOffset:  NSGAWiringBlinds.UnitSize
            );

            sim.RunWithSunVectorsBlinds(toSuns);
            sim.SaveLightHourGrid(outCsv);
        }
    }
}
