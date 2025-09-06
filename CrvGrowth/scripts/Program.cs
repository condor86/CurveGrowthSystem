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
        // —— 站点参数（与当前 LightingSimulator 默认一致：南京；如需更改请在此处改）——
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

            // === NSGA-II 基因边界（4 + 400 = 404）===
            const int repellerCount = 4;
            const int pointCount    = 400;
            const int geneLen       = repellerCount + pointCount;

            var lo = new double[geneLen];
            var hi = new double[geneLen];
            for (int i = 0; i < repellerCount; i++) { lo[i] = 0.01; hi[i] = 5.0; }
            for (int i = repellerCount; i < geneLen; i++) { lo[i] = 0.0; hi[i] = 100.0; }

            // === NSGA-II 日志目录（每代 front0 / bestGenes.csv）===
            string nsgaLogDir = Path.Combine(resultDir, "nsga_logs");
            Directory.CreateDirectory(nsgaLogDir);

            // === NSGA-II 配置：Evaluate 使用“预计算向量”的重载 ===
            var cfg = new NSGAConfig
            {
                PopulationSize = 50,
                Generations    = 100,
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

            // （可选）你的单步测试
            //TestSingleMoment.RunDefaultBatch();

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
            string outCrvCsv         = Path.Combine(resultDir, "resultsCrv.csv");
            string outLightingSummer = Path.Combine(resultDir, "resultsLighting_summer.csv");
            string outLightingWinter = Path.Combine(resultDir, "resultsLighting_winter.csv");

            SaveSolutionGeometryAndLighting(
                genes: rep.Genes,
                startingPoints: startingPoints,
                repellerPoints: repellerPoints,
                outCrvCsv: outCrvCsv,
                outLightingSummerCsv: outLightingSummer,
                outLightingWinterCsv: outLightingWinter,
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
            string outCrvCsv,
            string outLightingSummerCsv,
            string outLightingWinterCsv,
            Vector3[] summerToSuns,
            Vector3[] winterToSuns)
        {
            const int repellerCount = 4;
            const int offsetCount   = 400;

            // 1) 基因拆分
            var repellerFactors = genes.Take(repellerCount).ToList();
            var offsets         = genes.Skip(repellerCount).Take(offsetCount).ToArray();

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

            // 3) 转垂直（与你现有逻辑一致）：(x, y, 0) → (x, 0, z=y)
            var verticalCrv = flatCurve.Select(p => new Vector3(p.X, 0f, p.Y)).ToList();

            // 4) 逐点沿 -Y 偏移（前 N 个点；N = min(count, 400)）
            int N = Math.Min(verticalCrv.Count, offsets.Length);
            var extrudedCrv = new List<Vector3>(verticalCrv); // 独立列表，避免 alias
            for (int i = 0; i < N; i++)
            {
                var p = verticalCrv[i];
                extrudedCrv[i] = new Vector3(p.X, p.Y - (float)offsets[i], p.Z);
            }

            // 5) 导出竖直曲线（用于复盘/可视化）
            IOHelper.SavePointsToFile(outCrvCsv, extrudedCrv);

            // 6) 夏/冬分别用“向量直跑”并保存光照矩阵
            SimAndSaveVectors(verticalCrv, extrudedCrv, summerToSuns, outLightingSummerCsv);
            SimAndSaveVectors(verticalCrv, extrudedCrv, winterToSuns, outLightingWinterCsv);

            Console.WriteLine($"Saved: {outCrvCsv}");
            Console.WriteLine($"Saved: {outLightingSummerCsv}");
            Console.WriteLine($"Saved: {outLightingWinterCsv}");
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
                date:          NSGAWiring.SummerDate,  // 占位，不再用于太阳角计算
                startTime:     NSGAWiring.StartTime,
                endTime:       NSGAWiring.EndTime,
                interval:      NSGAWiring.Interval,
                roomWidth:     NSGAWiring.RoomWidth,
                roomDepth:     NSGAWiring.RoomDepth,
                gridSize:      NSGAWiring.GridSize
            );

            // 需要你在 LightingSimulator 中新增 RunWithSunVectors(Vector3[] toSuns) 方法
            sim.RunWithSunVectors(toSuns);
            sim.SaveLightHourGrid(outCsv);
        }
    }
}
