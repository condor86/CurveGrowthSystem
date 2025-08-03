using System;
using System.Collections.Generic;
using System.IO;
using System.Diagnostics;
using System.Linq;
using System.Numerics;

namespace CrvGrowth
{
    class Program
    {
        static void Main(string[] args)
        {
            var stopwatch1 = Stopwatch.StartNew();  // 开始计时
            
            // 根目录路径（即可执行文件所在目录）
            string rootDir = AppDomain.CurrentDomain.BaseDirectory;

            // 上级目录，用于保存输出结果
            string parentDir = Path.GetFullPath(Path.Combine(rootDir, "..", "..", ".."));

            // 输入文件路径（data/ 目录下）
            string startingPath = Path.Combine(parentDir, "data", "iStartingPositions.csv");
            string repellerPath = Path.Combine(parentDir, "data", "iRepellers.csv");
            string factorPath   = Path.Combine(parentDir, "data", "iRepellerFactors.csv");

            // 输出文件路径（results/ 目录下）
            string resultPathCrv      = Path.Combine(parentDir, "results", "resultsCrv.csv");
            string resultPathLighting = Path.Combine(parentDir, "results", "resultsLighting.csv");

            // 读取输入
            var startingPoints   = IOHelper.LoadPointsFromFile(startingPath);
            var repellerPoints   = IOHelper.LoadPointsFromFile(repellerPath);
            var repellerFactors  = IOHelper.LoadFactorsFromFile(factorPath);

            // 平面生长算法
            var system = new GrowthSystem();
            var flatCurve = system.Run(
                starting:        startingPoints,
                repellers:       repellerPoints,
                repellerFactors: repellerFactors,
                maxPointCount:   200,
                maxIterCount:    200,
                baseDist:        75.0
            );
            
            var resultCrv = flatCurve;

            // 保存平面曲线
            IOHelper.SavePointsToFile(resultPathCrv, resultCrv);
            Console.WriteLine($"共生成 {resultCrv.Count} 个点，结果已保存至：{resultPathCrv}");
            
            stopwatch1.Stop(); 
            Console.WriteLine($"Step1 平面生形耗时: {stopwatch1.ElapsedMilliseconds} ms");

            var stopwatch2 = Stopwatch.StartNew();

            // 生成完整几何体：将 XY 平面转为 XZ（Y → Z）
            List<Vector3> verticalCrv = flatCurve
                .Select(p => new Vector3(p.X, 0.0f, p.Y))
                .ToList();

            var extrudedCrv = new List<Vector3>(verticalCrv.Count);
            float offset = 100.0f;

            foreach (var pt in verticalCrv)
            {
                extrudedCrv.Add(new Vector3(pt.X, pt.Y - offset, pt.Z));
            }

            // 光照模拟
            var simulator = new LightingSimulator(
                verticalCurve: verticalCrv,
                extrudedCurve: extrudedCrv,
                date: new DateOnly(2025, 6, 25),
                startTime: new TimeOnly(8, 0),
                endTime: new TimeOnly(16, 0),
                interval: TimeSpan.FromHours(2),
                roomWidth: 1000.0,
                roomDepth: 1000.0,
                gridSize: 10.0
            );

            // 保存光照结果
            simulator.RunSimulation();
            simulator.SaveLightHourGrid(resultPathLighting);
            
            stopwatch2.Stop(); 
            Console.WriteLine($"Step2 光照模拟耗时: {stopwatch2.ElapsedMilliseconds} ms");
        }
    }
}
