// File: CrvGrowth/TestSingleMoment.cs
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Numerics;

namespace CrvGrowth
{
    /// <summary>
    /// 从 results/resultsCrv.csv 读取已“转竖直并挤出”的 extrudedCrv，
    /// 通过将 Y=0 的方式“反向挤出”得到 verticalCrv，
    /// 在给定时刻分别对夏/冬进行单时刻光照模拟，
    /// 将日照网格输出到 test_summer.csv / test_winter.csv。
    /// </summary>
    public static class TestSingleMoment
    {
        /// <param name="timeText">单次模拟的时刻（本地时区），格式 "HH:mm" ，如 "08:00"</param>
        /// <param name="inputExtrudedCsv">输入曲线（extrudedCrv）CSV；默认 results/resultsCrv.csv</param>
        /// <param name="outSummerCsv">夏季输出 CSV；默认 results/test_summer.csv</param>
        /// <param name="outWinterCsv">冬季输出 CSV；默认 results/test_winter.csv</param>
        public static void Run(
            string timeText = "08:00",
            string? inputExtrudedCsv = null,
            string? outSummerCsv = null,
            string? outWinterCsv = null
        )
        {
            // —— 目录与默认路径（与 Program.cs 一致）——
            string rootDir   = AppDomain.CurrentDomain.BaseDirectory;
            string parentDir = Path.GetFullPath(Path.Combine(rootDir, "..", "..", ".."));
            string resultDir = Path.Combine(parentDir, "results");
            Directory.CreateDirectory(resultDir);

            inputExtrudedCsv ??= Path.Combine(resultDir, "resultsCrv.csv");
            outSummerCsv     ??= Path.Combine(resultDir, "test_summer.csv");
            outWinterCsv     ??= Path.Combine(resultDir, "test_winter.csv");

            // —— 解析时刻 —— 
            if (!TimeSpan.TryParseExact(timeText, "hh\\:mm", CultureInfo.InvariantCulture, out var tod))
                throw new ArgumentException($"无法解析时刻：{timeText}（期望格式 HH:mm，例如 08:00）");
            var tOnly = new TimeOnly(tod.Hours, tod.Minutes);

            // —— 读取 extrudedCrv —— 
            if (!File.Exists(inputExtrudedCsv))
                throw new FileNotFoundException($"未找到输入曲线文件：{inputExtrudedCsv}");
            var extrudedCrv = LoadPointsFlexible(inputExtrudedCsv);
            if (extrudedCrv.Count == 0)
                throw new InvalidOperationException($"输入曲线为空：{inputExtrudedCsv}");

            // —— 反向挤出：得到 verticalCrv（投影到 XZ 平面：Y=0）——
            var verticalCrv = new List<Vector3>(extrudedCrv.Count);
            for (int i = 0; i < extrudedCrv.Count; i++)
            {
                var p = extrudedCrv[i];
                verticalCrv.Add(new Vector3(p.X, 0f, p.Z));
            }

            // —— 夏至 / 冬至：单时刻模拟并保存 —— 
            SimOneMoment(verticalCrv, extrudedCrv, NSGAWiring.SummerDate, tOnly, outSummerCsv);
            SimOneMoment(verticalCrv, extrudedCrv, NSGAWiring.WinterDate, tOnly, outWinterCsv);

            Console.WriteLine($"[TestSingleMoment] Done. time={timeText}");
            Console.WriteLine($"  Input(extruded): {inputExtrudedCsv}");
            Console.WriteLine($"  Summer grid  ->  {outSummerCsv}");
            Console.WriteLine($"  Winter grid  ->  {outWinterCsv}");
        }

        /// <summary>
        /// 单时刻模拟：start==end==tOnly，interval=1min。
        /// 其它参数复用 NSGAWiring 的全局设定。
        /// </summary>
        private static void SimOneMoment(
            List<Vector3> verticalCrv,
            List<Vector3> extrudedCrv,
            DateOnly date,
            TimeOnly tOnly,
            string outCsv)
        {
            var sim = new LightingSimulator(
                verticalCurve: verticalCrv,
                extrudedCurve: extrudedCrv,
                date:          date,
                startTime:     tOnly,
                endTime:       tOnly,
                interval:      TimeSpan.FromHours(2),
                roomWidth:     NSGAWiring.RoomWidth,
                roomDepth:     NSGAWiring.RoomDepth,
                gridSize:      NSGAWiring.GridSize
            );

            sim.RunSimulation();
            sim.SaveLightHourGrid(outCsv);
            Console.WriteLine($"Saved: {outCsv}");
        }

        /// <summary>
        /// 兼容两种常见写法：
        /// 1) 每行一个点：x,y,z
        /// 2) 一行多个三元组：x,y,z x,y,z ...
        /// 允许逗号与空白混合分隔。
        /// </summary>
        private static List<Vector3> LoadPointsFlexible(string path)
        {
            var text = File.ReadAllText(path);
            var seps = new[] { ' ', '\t', '\r', '\n', ',' };
            var toks = text.Split(seps, StringSplitOptions.RemoveEmptyEntries);

            var pts = new List<Vector3>(Math.Max(16, toks.Length / 3));
            var n = toks.Length - (toks.Length % 3);

            for (int i = 0; i < n; i += 3)
            {
                if (double.TryParse(toks[i],   NumberStyles.Float, CultureInfo.InvariantCulture, out var x) &&
                    double.TryParse(toks[i+1], NumberStyles.Float, CultureInfo.InvariantCulture, out var y) &&
                    double.TryParse(toks[i+2], NumberStyles.Float, CultureInfo.InvariantCulture, out var z))
                {
                    pts.Add(new Vector3((float)x, (float)y, (float)z));
                }
            }
            return pts;
        }
    }
}
