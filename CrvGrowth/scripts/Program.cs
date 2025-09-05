// File: CrvGrowth/Program.cs
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Numerics;
using NSGAII;   // 需要你已经添加 NSGAII.cs
using CrvGrowth; // 与当前命名空间一致

namespace CrvGrowth
{
    class Program
    {
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

            // === NSGA-II 基因边界（4 + 400 = 404）===
            const int repellerCount = 4;
            const int pointCount    = 400;
            const int geneLen       = repellerCount + pointCount;

            var lo = new double[geneLen];
            var hi = new double[geneLen];
            for (int i = 0; i < repellerCount; i++) { lo[i] = 0.01; hi[i] = 5.0;   }   // 4 个 repeller 因子
            for (int i = repellerCount; i < geneLen; i++) { lo[i] = 0.0;  hi[i] = 100.0; } // 400 个逐点位移（沿 -Y）

            // === NSGA-II 日志目录（每代 front0 / bestGenes.csv）===
            string nsgaLogDir = Path.Combine(resultDir, "nsga_logs");
            Directory.CreateDirectory(nsgaLogDir);

            // === NSGA-II 配置===
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

                // 使用 NSGAWiring 连接 GrowthSystem 与 LightingSimulator
                Evaluate = NSGAWiring.MakeEvaluator(startingPoints, repellerPoints)
                // 如需从文件懒加载（首轮读取，后续复用）：
                // Evaluate = NSGAWiring.MakeEvaluatorFromFiles(startingCsv, repellerCsv)
            };
            
            
            TestSingleMoment.RunDefaultBatch();
            
            
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

            // === 导出代表解的几何与光照（夏/冬各一份）===
            string outCrvCsv         = Path.Combine(resultDir, "resultsCrv.csv");
            string outLightingSummer = Path.Combine(resultDir, "resultsLighting_summer.csv");
            string outLightingWinter = Path.Combine(resultDir, "resultsLighting_winter.csv");

                SaveSolutionGeometryAndLighting(
                genes: rep.Genes,
                startingPoints: startingPoints,
                repellerPoints: repellerPoints,
                outCrvCsv: outCrvCsv,
                outLightingSummerCsv: outLightingSummer,
                outLightingWinterCsv: outLightingWinter
            );

            totalWatch.Stop();
            Console.WriteLine($"All done. Total time: {totalWatch.Elapsed}");
        }

        /// 导出解的函数
        private static void SaveSolutionGeometryAndLighting(
            double[] genes,
            List<Vector3> startingPoints,
            List<Vector3> repellerPoints,
            string outCrvCsv,
            string outLightingSummerCsv,
            string outLightingWinterCsv)
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

            // 6) 导出竖直曲线（用于复盘/可视化）
            IOHelper.SavePointsToFile(outCrvCsv, extrudedCrv);

            // 7) 夏/冬分别跑一次并保存光照矩阵
            SimAndSave(verticalCrv, extrudedCrv, NSGAWiring.SummerDate, outLightingSummerCsv);
            SimAndSave(verticalCrv, extrudedCrv, NSGAWiring.WinterDate, outLightingWinterCsv);

            Console.WriteLine($"Saved: {outCrvCsv}");
            Console.WriteLine($"Saved: {outLightingSummerCsv}");
            Console.WriteLine($"Saved: {outLightingWinterCsv}");
        }

        private static void SimAndSave(
            List<Vector3> verticalCrv,
            List<Vector3> extrudedCrv,
            DateOnly date,
            string outCsv)
        {
            var sim = new LightingSimulator(
                verticalCurve: verticalCrv,
                extrudedCurve: extrudedCrv,
                date:          date,
                startTime:     NSGAWiring.StartTime,
                endTime:       NSGAWiring.EndTime,
                interval:      NSGAWiring.Interval,
                roomWidth:     NSGAWiring.RoomWidth,
                roomDepth:     NSGAWiring.RoomDepth,
                gridSize:      NSGAWiring.GridSize
            );
            sim.RunSimulation();
            sim.SaveLightHourGrid(outCsv);
        }
    }
}
