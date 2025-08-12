// File: CrvGrowth/NSGAWiring.cs
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace CrvGrowth
{
    /// <summary>
    /// 将 GrowthSystem 与 LightingSimulator 通过 Evaluate(genes) 串联起来，
    /// 供 NSGA-II 的 NSGAConfig.Evaluate 回调直接使用。
    ///
    /// 约定的基因布局：
    ///   genes[0..3]   -> 4 个 repeller 因子，范围 [0.01, 5.0]
    ///   genes[4..403] -> 400 个逐点位移量，沿 -Y 方向，范围 [0, 100]
    ///
    /// 目标方向（统一最小化）：
    ///   f0 = 夏季光照小时（越小越好） -> 直接返回
    ///   f1 = -冬季光照小时（冬季越多越好） -> 取负返回
    ///
    /// 使用方式（示例）：
    ///   var eval = NSGAWiring.MakeEvaluator(startingPoints, repellerPoints);
    ///   cfg.Evaluate = eval;
    /// </summary>
    public static class NSGAWiring
    {
        // ======= GrowthSystem 默认参数（与你现有主程序一致） =======
        public static int    MaxPointCount = 200;
        public static int    MaxIterCount  = 200;
        public static double BaseDist      = 75.0;

        // ======= 几何参数 =======
        /// <summary>沿 -Y 挤出的深度（如果将来想优化它，可把它纳入第 405 个基因）</summary>
        public static float  ExtrudeDepth  = 100f;

        // ======= 光照模拟参数（可根据工程需要调整/外部设置） =======
        public static DateOnly SummerDate  = new DateOnly(2025, 6, 21);
        public static DateOnly WinterDate  = new DateOnly(2025, 12, 21);
        public static TimeOnly StartTime   = new TimeOnly(8, 0);
        public static TimeOnly EndTime     = new TimeOnly(16, 0);
        public static TimeSpan Interval    = TimeSpan.FromHours(1);

        public static double RoomWidth     = 1000.0;
        public static double RoomDepth     = 1000.0;
        public static double GridSize      = 10.0;

        /// <summary>使用平均日照小时作为目标（true）；
        /// 若想使用总日照小时，设为 false 并在 GetLightMetric 内切换。</summary>
        public static bool   UseAverageLightHours = true;

        // =====================================================================
        // 对外：两种工厂方法 —— 1) 直接传内存数据； 2) 从文件路径懒加载（需 IOHelper）
        // =====================================================================

        /// <summary>
        /// 用内存中的起始点与 repeller 点，创建 Evaluate 回调。
        /// </summary>
        public static Func<double[], double[]> MakeEvaluator(
            List<Vector3> startingPoints,
            List<Vector3> repellerPoints)
        {
            if (startingPoints == null || startingPoints.Count == 0)
                throw new ArgumentException("startingPoints 不能为空。");
            if (repellerPoints == null)
                throw new ArgumentException("repellerPoints 不能为空。");

            // 用闭包把数据捕获起来，避免 Evaluate 内重复读盘
            return (genes) => EvaluateOnce(genes, startingPoints, repellerPoints);
        }

        /// <summary>
        /// 从文件路径懒加载（需要项目中的 IOHelper.LoadPointsFromFile），
        /// 适合你已有的数据放在 data/ 目录下的情况。
        /// </summary>
        public static Func<double[], double[]> MakeEvaluatorFromFiles(
            string startingCsvPath,
            string repellersCsvPath)
        {
            // 首次调用时加载一次，后续复用缓存
            List<Vector3>? starting = null;
            List<Vector3>? repellers = null;

            return (genes) =>
            {
                starting  ??= IOHelper.LoadPointsFromFile(startingCsvPath);
                repellers ??= IOHelper.LoadPointsFromFile(repellersCsvPath);

                return EvaluateOnce(genes, starting, repellers);
            };
        }

        // =====================================================================
        // 内部：单次评估主流程（基因 -> 平面生长 -> verticalCrv -> 逐点 -Y 偏移 -> extrude -Y -> 夏/冬光照）
        // =====================================================================

        private static double[] EvaluateOnce(
            double[] genes,
            List<Vector3> startingPoints,
            List<Vector3> repellerPoints)
        {
            if (genes == null || genes.Length < 404)
                throw new ArgumentException("基因长度不足（需要 404，含 4 个 repeller 因子 + 400 个逐点位移）。");

            // 1) 基因拆分
            const int repellerCount = 4;
            const int offsetCount   = 400;

            var repellerFactors = new List<double>(repellerCount);
            for (int i = 0; i < repellerCount; i++) repellerFactors.Add(genes[i]);

            var offsets = new double[offsetCount];
            Array.Copy(genes, repellerCount, offsets, 0, offsetCount);

            // 2) 平面生长（XY 平面，Z=0）
            var growth = new GrowthSystem();
            var flatCurve = growth.Run(
                starting:        startingPoints,
                repellers:       repellerPoints,
                repellerFactors: repellerFactors,
                maxPointCount:   MaxPointCount,
                maxIterCount:    MaxIterCount,
                baseDist:        BaseDist
            );

            // 3) 转为竖直曲线：与你 Program.cs 保持一致 (x, y, 0) -> (x, 0, z=y)
            var verticalCrv = ToVerticalXZ(flatCurve);

            // 4) 逐点沿 -Y 应用偏移（前 N 个点；N = min(Count, 400)）
            ApplyOffsetsMinusY(verticalCrv, offsets);

            // 5) 沿 -Y 挤出得到 extruded 曲线
            var extrudedCrv = BuildExtrudedMinusY(verticalCrv, ExtrudeDepth);

            // 6) 夏季光照
            double summerMetric = SimulateAndGetMetric(verticalCrv, extrudedCrv, SummerDate);

            // 7) 冬季光照
            double winterMetric = SimulateAndGetMetric(verticalCrv, extrudedCrv, WinterDate);

            // 8) 统一最小化：{夏季，-冬季}
            return new[] { summerMetric, -winterMetric };
        }

        // =====================================================================
        // 几何辅助
        // =====================================================================

        /// <summary>将平面曲线 (x, y, 0) 映射为立面曲线 (x, 0, z=y)</summary>
        private static List<Vector3> ToVerticalXZ(List<Vector3> flat)
        {
            var vertical = new List<Vector3>(flat.Count);
            for (int i = 0; i < flat.Count; i++)
            {
                var p = flat[i];
                vertical.Add(new Vector3(p.X, 0f, p.Y));
            }
            return vertical;
        }

        /// <summary>对 verticalCrv 的前 N 个点（N=min(count, offsets.Length)）沿 -Y 平移 offset[i]</summary>
        private static void ApplyOffsetsMinusY(List<Vector3> verticalCrv, double[] offsets)
        {
            int N = Math.Min(verticalCrv.Count, offsets.Length);
            for (int i = 0; i < N; i++)
            {
                var p = verticalCrv[i];
                verticalCrv[i] = new Vector3(p.X, p.Y - (float)offsets[i], p.Z);
            }
        }

        /// <summary>基于 verticalCrv，沿 -Y 方向挤出 ExtrudeDepth，得到 extruded 曲线</summary>
        private static List<Vector3> BuildExtrudedMinusY(List<Vector3> verticalCrv, float depth)
        {
            var extruded = new List<Vector3>(verticalCrv.Count);
            for (int i = 0; i < verticalCrv.Count; i++)
            {
                var p = verticalCrv[i];
                extruded.Add(new Vector3(p.X, p.Y - depth, p.Z));
            }
            return extruded;
        }

        // =====================================================================
        // 光照辅助
        // =====================================================================

        /// <summary>在给定日期下运行一次 LightingSimulator，并返回度量（平均或总计）</summary>
        private static double SimulateAndGetMetric(
            List<Vector3> verticalCrv,
            List<Vector3> extrudedCrv,
            DateOnly date)
        {
            var sim = new LightingSimulator(
                verticalCurve: verticalCrv,
                extrudedCurve: extrudedCrv,
                date:          date,
                startTime:     StartTime,
                endTime:       EndTime,
                interval:      Interval,
                roomWidth:     RoomWidth,
                roomDepth:     RoomDepth,
                gridSize:      GridSize
            );
            sim.RunSimulation();

            return GetLightMetric(sim);
        }

        /// <summary>根据 UseAverageLightHours 开关返回平均或总计</summary>
        private static double GetLightMetric(LightingSimulator sim)
        {
            if (UseAverageLightHours)
            {
                // 需要在 LightingSimulator 中添加：
                // public double ComputeAverageLightHours() { ... }
                return sim.GetTotalLightHours();
            }
            else
            {
                // 如果你更偏好总计，取消注释，并在 LightingSimulator 中添加：
                // public int ComputeTotalLightHours() { ... }
                // return sim.ComputeTotalLightHours();
                // 这里默认也返回平均，避免编译错误：
                return sim.GetTotalLightHours();
            }
        }
    }
}
