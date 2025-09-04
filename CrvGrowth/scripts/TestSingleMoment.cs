// File: CrvGrowth/TestSingleMoment.cs
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Numerics;

namespace CrvGrowth
{
    public static class TestSingleMoment
    {
        /// <summary>
        /// 单时刻：从 results/resultsCrv.csv 读取 extrudedCrv（已转竖直且挤出），
        /// 反向挤出得到 verticalCrv（Y=0 → 投影到 XZ），
        /// 在给定时刻分别对夏至/冬至进行单时刻光照模拟，
        /// 输出到 {prefix}_{h[_mm]}_summer.csv / {prefix}_{h[_mm]}_winter.csv。
        /// </summary>
        public static void Run(
            string timeText = "08:00",
            string? inputExtrudedCsv = null,
            string? outDir = null,
            string prefix = "test"
        )
        {
            // 路径与 Program.cs 保持一致
            string rootDir   = AppDomain.CurrentDomain.BaseDirectory;
            string parentDir = Path.GetFullPath(Path.Combine(rootDir, "..", "..", ".."));
            string resultDir = Path.Combine(parentDir, "results");

            inputExtrudedCsv ??= Path.Combine(resultDir, "resultsCrv.csv");
            outDir           ??= Path.Combine(resultDir, "tests_single_moment");
            Directory.CreateDirectory(outDir);

            // 解析时刻
            if (!TimeSpan.TryParseExact(timeText, "hh\\:mm", CultureInfo.InvariantCulture, out var tod))
                throw new ArgumentException($"无法解析时刻：{timeText}（期望格式 HH:mm，例如 08:00）");
            var tOnly  = new TimeOnly(tod.Hours, tod.Minutes);
            var label  = FormatTimeLabel(tOnly);
            var outSum = Path.Combine(outDir, $"{prefix}_{label}_summer.csv");
            var outWin = Path.Combine(outDir, $"{prefix}_{label}_winter.csv");

            // 读取 extrudedCrv（兼容一行多个三元组或逐行三元组）
            if (!File.Exists(inputExtrudedCsv))
                throw new FileNotFoundException($"未找到输入曲线文件：{inputExtrudedCsv}");
            var extrudedCrv = LoadPointsFlexible(inputExtrudedCsv);
            if (extrudedCrv.Count == 0)
                throw new InvalidOperationException($"输入曲线为空：{inputExtrudedCsv}");

            // 反向挤出：verticalCrv = (x, 0, z)
            var verticalCrv = new List<Vector3>(extrudedCrv.Count);
            for (int i = 0; i < extrudedCrv.Count; i++)
            {
                var p = extrudedCrv[i];
                verticalCrv.Add(new Vector3(p.X, 0f, p.Z));
            }

            // 夏至 / 冬至：单时刻模拟并保存（start==end==tOnly，interval=1min）
            SimOneMoment(verticalCrv, extrudedCrv, NSGAWiring.SummerDate, tOnly, outSum);
            SimOneMoment(verticalCrv, extrudedCrv, NSGAWiring.WinterDate, tOnly, outWin);

            Console.WriteLine($"[TestSingleMoment] Done. time={timeText}");
            Console.WriteLine($"  Input(extruded): {inputExtrudedCsv}");
            Console.WriteLine($"  Summer grid  ->  {outSum}");
            Console.WriteLine($"  Winter grid  ->  {outWin}");
        }

        /// <summary>
        /// 批量：对一组时刻分别输出夏/冬结果。
        /// 仅加载一次 extrudedCrv 并生成一次 verticalCrv，随后循环模拟以减少 I/O。
        /// </summary>
        public static void RunBatch(
            IEnumerable<TimeOnly> times,
            string? inputExtrudedCsv = null,
            string? outDir = null,
            string prefix = "test"
        )
        {
            // 路径与 Program.cs 保持一致
            string rootDir   = AppDomain.CurrentDomain.BaseDirectory;
            string parentDir = Path.GetFullPath(Path.Combine(rootDir, "..", "..", ".."));
            string resultDir = Path.Combine(parentDir, "results");

            inputExtrudedCsv ??= Path.Combine(resultDir, "resultsCrv.csv");
            outDir           ??= Path.Combine(resultDir, "tests_single_moment");
            Directory.CreateDirectory(outDir);

            if (!File.Exists(inputExtrudedCsv))
                throw new FileNotFoundException($"未找到输入曲线文件：{inputExtrudedCsv}");
            var extrudedCrv = LoadPointsFlexible(inputExtrudedCsv);
            if (extrudedCrv.Count == 0)
                throw new InvalidOperationException($"输入曲线为空：{inputExtrudedCsv}");

            // 反向挤出：verticalCrv = (x, 0, z)
            var verticalCrv = new List<Vector3>(extrudedCrv.Count);
            for (int i = 0; i < extrudedCrv.Count; i++)
            {
                var p = extrudedCrv[i];
                verticalCrv.Add(new Vector3(p.X, 0f, p.Z));
            }

            foreach (var t in times)
            {
                var label  = FormatTimeLabel(t);
                var outSum = Path.Combine(outDir, $"{prefix}_{label}_summer.csv");
                var outWin = Path.Combine(outDir, $"{prefix}_{label}_winter.csv");

                SimOneMoment(verticalCrv, extrudedCrv, NSGAWiring.SummerDate, t, outSum);
                SimOneMoment(verticalCrv, extrudedCrv, NSGAWiring.WinterDate, t, outWin);

                Console.WriteLine($"[TestSingleMoment] time={t} -> {Path.GetFileName(outSum)}, {Path.GetFileName(outWin)}");
            }
        }

        /// <summary>
        /// 默认批量：08:00, 10:00, 12:00, 14:00, 16:00。
        /// </summary>
        public static void RunDefaultBatch(
            string? inputExtrudedCsv = null,
            string? outDir = null,
            string prefix = "test"
        )
        {
            var times = new[]
            {
                new TimeOnly(8, 0),
                new TimeOnly(10, 0),
                new TimeOnly(12, 0),
                new TimeOnly(14, 0),
                new TimeOnly(16, 0),
            };
            RunBatch(times, inputExtrudedCsv, outDir, prefix);
        }

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
                interval:      TimeSpan.FromMinutes(1),
                roomWidth:     NSGAWiring.RoomWidth,
                roomDepth:     NSGAWiring.RoomDepth,
                gridSize:      NSGAWiring.GridSize
            );
            sim.RunSimulation();
            sim.SaveLightHourGrid(outCsv);
        }

        // 文件读取：支持逐行 x,y,z 或一行多个 x,y,z 三元组（逗号/空白混合分隔）
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

        // 文件名中的时刻标签：8:00 -> "8"，8:30 -> "8_30"
        private static string FormatTimeLabel(TimeOnly t)
            => t.Minute == 0 ? $"{t.Hour}" : $"{t.Hour}_{t.Minute}";
    }
}
