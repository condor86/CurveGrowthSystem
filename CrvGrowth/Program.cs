using System;
using System.Collections.Generic;
using System.IO;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using NumSharp;

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

            // 输入文件路径
            string startingPath  = Path.Combine(rootDir, "iStartingPositions.txt");
            string repellerPath  = Path.Combine(rootDir, "iRepellers.txt");
            string factorPath    = Path.Combine(rootDir, "iRepellerFactors.txt");

            // 输出文件路径
            string resultPathCrv      = Path.Combine(parentDir, "resultsCrv.txt");
            string resultPathLighting = Path.Combine(parentDir, "resultsLighting.txt");

            // ==== 读取输入 ====
            NDArray startingND  = IOHelper.LoadPointsAsNDArray(startingPath);
            NDArray repellerND  = IOHelper.LoadPointsAsNDArray(repellerPath);
            NDArray factorND    = IOHelper.LoadFactorsAsNDArray(factorPath);

            // ==== 执行平面生长 ====
            var system = new GrowthSystem();
            NDArray resultND = system.Run(
                starting:        startingND,
                repellers:       repellerND,
                repellerFactors: factorND,
                maxPointCount:   200,
                maxIterCount:    200,
                baseDist:        75.0
            );

            // ==== 保存结果 ====
            IOHelper.SaveNDArrayAsPointFile(resultPathCrv, resultND);
            Console.WriteLine($"共生成 {resultND.shape[0]} 个点，结果已保存至：{resultPathCrv}");

            stopwatch1.Stop(); 
            Console.WriteLine($"Step1 平面生形耗时: {stopwatch1.ElapsedMilliseconds} ms");

            // ==== 光照模拟部分 ====
            var stopwatch2 = Stopwatch.StartNew();

            // 将 XY 平面点转换为 XZ 垂直曲线点（中间用 Vector3 桥接）
            var verticalCrv = new List<System.Numerics.Vector3>();
            for (int i = 0; i < resultND.shape[0]; i++)
            {
                float x = (float)resultND[i, 0];
                float z = (float)resultND[i, 1];  // 注意：Y 轴 → Z
                verticalCrv.Add(new System.Numerics.Vector3(x, 0.0f, z));
            }

            // 向下 extrude 一段高度
            var extrudedCrv = new List<System.Numerics.Vector3>();
            float offset = 100.0f;
            foreach (var pt in verticalCrv)
            {
                extrudedCrv.Add(new System.Numerics.Vector3(pt.X, pt.Y - offset, pt.Z));
            }

            // ==== 光照模拟 ====
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

            simulator.RunSimulation();
            simulator.SaveLightHourGrid(resultPathLighting);

            stopwatch2.Stop(); 
            Console.WriteLine($"Step2 光照模拟耗时: {stopwatch2.ElapsedMilliseconds} ms");
        }
    }
}
