// File: CrvGrowth/Program.cs
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Numerics;

using NSGAII;                // 需要 NSGAII.cs
using CrvGrowth;             // 当前命名空间
using CrvGrowth.Scripts;     // scripts/SunCache.cs 里的 SunVectors

namespace CrvGrowth
{
    class Program
    {
        // —— 站点参数（与 LightingSimulator/NSGAWiring 默认一致：南京；如需更改请在此处改）——
        private const double SiteLatitudeDeg   = 32.0603;
        private const double SiteLongitudeDeg  = 118.7969;
        private const double SiteTimezoneHours = 8.0;

        // —— 模型坐标系：Up=+Z, North=+Y（如有模型相对真北的偏航，可在此处旋转 North）——
        private static readonly Vector3 Up    = new(0, 0, 1);
        private static readonly Vector3 North = new(0, 1, 0);

        static void Main(string[] args)
        {
            var totalWatch = Stopwatch.StartNew();

            // === 目录结构 ===
            string rootDir   = AppDomain.CurrentDomain.BaseDirectory;
            string parentDir = Path.GetFullPath(Path.Combine(rootDir, "..", "..", ".."));
            string dataDir   = Path.Combine(parentDir, "data");
            string resultDir = Path.Combine(parentDir, "results");
            Directory.CreateDirectory(resultDir);

            // === 输入数据路径 ===
            string startingCsv = Path.Combine(dataDir, "iStartingPositions.csv");
            string repellerCsv = Path.Combine(dataDir, "iRepellers.csv");

            // === 预加载输入（避免 Evaluate 内反复读盘）===
            var startingPoints = IOHelper.LoadPointsFromFile(startingCsv);
            var repellerPoints = IOHelper.LoadPointsFromFile(repellerCsv);

            // === 预计算太阳向量（一次）：夏/冬各一组，包含 [StartTime, EndTime] 的全部采样点 ===
            var summerToSuns = SunVectors.Build(
                NSGAWiring.SummerDate, NSGAWiring.StartTime, NSGAWiring.EndTime, NSGAWiring.Interval,
                SiteLatitudeDeg, SiteLongitudeDeg, SiteTimezoneHours,
                Up, North);

            var winterToSuns = SunVectors.Build(
                NSGAWiring.WinterDate, NSGAWiring.StartTime, NSGAWiring.EndTime, NSGAWiring.Interval,
                SiteLatitudeDeg, SiteLongitudeDeg, SiteTimezoneHours,
                Up, North);

            // === NSGA-II 基因边界（4 + 400 + 400 = 804）===
            const int repellerCount = 4;     // 斥力因子
            const int offsetCount   = 400;   // 逐点 -Y 偏移
            const int angleCount    = 400;   // 逐点局部平面旋转角
            const int geneLen       = repellerCount + offsetCount + angleCount;

            var lo = new double[geneLen];
            var hi = new double[geneLen];

            // 1) 4 个 repeller 因子
            for (int i = 0; i < repellerCount; i++) { lo[i] = 0.01; hi[i] = 5.0; }

            // 2) 400 个逐点 -Y 偏移（与你当前设置保持一致：50~100）
            for (int i = repellerCount; i < repellerCount + offsetCount; i++) { lo[i] = 50.0; hi[i] = 100.0; }

            // 3) 400 个逐点旋转角（单位：度）
            for (int i = repellerCount + offsetCount; i < geneLen; i++) { lo[i] = -25.0; hi[i] = 25.0; }

            // === NSGA-II 日志目录（每代 front0 / bestGenes.csv）===
            string nsgaLogDir = Path.Combine(resultDir, "nsga_logs");
            Directory.CreateDirectory(nsgaLogDir);

            // === NSGA-II 配置：Evaluate 使用“预计算向量”的重载 ===
            var cfg = new NSGAConfig
            {
                PopulationSize = 50,
                Generations    = 1,
                CrossoverRate  = 0.9,
                MutationRate   = 1.0 / geneLen,  // ≈ 1/n
                GeneLength     = geneLen,
                LowerBounds    = lo,
                UpperBounds    = hi,
                RandomSeed     = 1,
                DegreeOfParallelism = Environment.ProcessorCount,
                LogDir         = nsgaLogDir,
                SbxEta         = 20.0,
                PolyMutationEta= 20.0,

                // 关键：把夏/冬两组“指向太阳”的单位向量数组传给 NSGAWiring
                Evaluate = NSGAWiring.MakeEvaluator(
                    startingPoints, repellerPoints,
                    summerToSuns, winterToSuns)
            };

            // === 运行 NSGA-II ===
            Console.WriteLine("NSGA-II optimization started...");
            var runWatch = Stopwatch.StartNew();

            var solver   = new NSGAII.NSGAII(cfg);
            var finalPop = solver.Run();

            runWatch.Stop();
            Console.WriteLine($"NSGA-II finished in {runWatch.Elapsed}.");

            // === 取最终一代 Pareto 前沿 & 代表解（按目标向量 L1 和）===
            var pareto = finalPop.Where(ind => ind.Rank == 0).ToList();
            Console.WriteLine($"Final Pareto size: {pareto.Count}");

            var rep = pareto.OrderBy(ind => ind.Objectives.Sum()).First();
            Console.WriteLine("Exporting representative solution geometry & lighting...");

            // === 导出代表解的几何与光照（同样使用“预计算向量”）===
            string outCrvCsv              = Path.Combine(resultDir, "resultsCrv.csv");
            string outVerticalCsv         = Path.Combine(resultDir, "resultsVertical.csv");          // verticalCrv
            string outFilletVerticalCsv   = Path.Combine(resultDir, "resultsFilletVertical.csv");    // verticalCrv 圆角化
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

            totalWatch.Stop();
            Console.WriteLine($"All done. Total time: {totalWatch.Elapsed}");
        }

        /// 导出解：几何与光照（光照使用“向量直跑”，保持与优化一致）
        private static void SaveSolutionGeometryAndLighting(
            double[] genes,
            List<Vector3> startingPoints,
            List<Vector3> repellerPoints,
            string outVerticalCsv,               // verticalCrv
            string outFilletVerticalCsv,         // verticalCrv 的圆角化版本
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

            // 1) 基因拆分
            var repellerFactors = genes.Take(repellerCount).ToList();
            var offsets         = genes.Skip(repellerCount).Take(offsetCount).ToArray();
            var anglesDeg       = genes.Skip(repellerCount + offsetCount).Take(angleCount).ToArray();

            // 2) 平面生长
            var system = new GrowthSystem();
            var flatCurve = system.Run(
                starting:        startingPoints,
                repellers:       repellerPoints,
                repellerFactors: repellerFactors,
                maxPointCount:   NSGAWiring.MaxPointCount,
                maxIterCount:    NSGAWiring.MaxIterCount,
                baseDist:        NSGAWiring.BaseDist
            );

            // 3) 转垂直（(x, y, 0) → (x, 0, z=y)）
            var verticalCrv = flatCurve.Select(p => new Vector3(p.X, 0f, p.Y)).ToList();

            // —— 导出 verticalCrv —— 
            IOHelper.SavePointsToFile(outVerticalCsv, verticalCrv);

            // —— 导出 verticalCrv 的“圆角化版本”（固定 9 采样点；是否闭合可按需要调整）——
            var filletVertical = FilletUtil.FilletPolylineWithFixedArcPoints(
                pts: verticalCrv,
                radius: filletRadiusDefault,
                arcPointCount: 9,
                angleEpsDeg: 1.0f,
                isClosed: true,          // 如 verticalCrv 是开口折线，将其改为 false
                clampRadius: true);
            IOHelper.SavePointsToFile(outFilletVerticalCsv, filletVertical);

            // 4) 逐点沿 -Y 偏移（前 N 个点）
            int N = Math.Min(verticalCrv.Count, offsets.Length);
            var extrudedCrv = new List<Vector3>(verticalCrv); // 独立列表，避免 alias
            for (int i = 0; i < N; i++)
            {
                var p = verticalCrv[i];
                extrudedCrv[i] = new Vector3(p.X, p.Y - (float)offsets[i], p.Z);
            }

            // 5) 逐点“局部平面旋转”（以 Pn 为枢轴，在由角平分线与 pnPn 线生成的平面内）
            NSGAWiring.ApplyLocalPlaneRotation(verticalCrv, extrudedCrv, anglesDeg);

            // 6) 导出挤出（含旋转）后的曲线（用于复盘/可视化）
            IOHelper.SavePointsToFile(outCrvCsv, extrudedCrv);

            // 7) 夏/冬分别用“向量直跑”并保存光照矩阵（传入镜像模式）
            SimAndSaveVectors(verticalCrv, extrudedCrv, summerToSuns, outLightingSummerCsv);
            SimAndSaveVectors(verticalCrv, extrudedCrv, winterToSuns, outLightingWinterCsv);

            Console.WriteLine($"Saved: {outVerticalCsv}");
            Console.WriteLine($"Saved: {outFilletVerticalCsv}");
            Console.WriteLine($"Saved: {outCrvCsv}");
            Console.WriteLine($"Saved: {outLightingSummerCsv}");
            Console.WriteLine($"Saved: {outLightingWinterCsv}");
            
            /*
            // 可选：NURBS 采样导出
            var extrudedNurbsCurve = NurbsTools.BuildRhinoLikeCurve(extrudedCrv, degree: 3);
            var ePts400 = NurbsTools.SampleByArcLength(extrudedNurbsCurve, count: 400);
            IOHelper.SavePointsToFile(outNurbsCsv, ePts400);
            Console.WriteLine($"Saved: {outNurbsCsv}");
            */
            
            var filletSampled = FilletUtil.FilletPolylineWithFixedArcPoints(
                pts: extrudedCrv,
                radius: filletRadiusDefault,
                arcPointCount: 9,      // 固定 9 点
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
                date:          NSGAWiring.SummerDate,  // 占位，不在“向量直跑”中使用
                startTime:     NSGAWiring.StartTime,
                endTime:       NSGAWiring.EndTime,
                interval:      NSGAWiring.Interval,
                roomWidth:     NSGAWiring.RoomWidth,
                roomDepth:     NSGAWiring.RoomDepth,
                gridSize:      NSGAWiring.GridSize,
                isClosed:      true,
                mirrorMode:    NSGAWiring.MirrorMode   // 关键：统一镜像模式
            );

            sim.RunWithSunVectors(toSuns);
            sim.SaveLightHourGrid(outCsv);
        }
    }
}
