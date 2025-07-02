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

            // 输入文件路径（位于可执行目录中）
            string startingPath  = Path.Combine(rootDir, "iStartingPositions.txt");
            string repellerPath  = Path.Combine(rootDir, "iRepellers.txt");
            string factorPath    = Path.Combine(rootDir, "iRepellerFactors.txt");

            // 输出文件路径（保存到工程目录根部）
            string resultPathCrv      = Path.Combine(parentDir, "resultsCrv.txt");
            string resultPathLighting = Path.Combine(parentDir, "resultsLighting.txt");

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
            Console.WriteLine($"共生成 {resultCrv.Count()} 个点，结果已保存至：{resultPathCrv}");
            
            stopwatch1.Stop(); 
            Console.WriteLine($"Step1 平面生形耗时: {stopwatch1.ElapsedMilliseconds} ms");

            var stopwatch2 = Stopwatch.StartNew();

            // 生成完整几何体
            List<Point3D> verticalCrv = flatCurve.Select(p => p.RotateXYToXZ()).ToList();

            var extrudedCrv = new List<Point3D>(verticalCrv.Count);
            double offset = 100.0;
            
            foreach (var pt in verticalCrv)
            {
                extrudedCrv.Add(new Point3D(pt.X, pt.Y - offset, pt.Z));
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
