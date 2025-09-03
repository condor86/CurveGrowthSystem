// File: CrvGrowth/NSGAWiring.cs
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Threading;

namespace CrvGrowth
{
    /// <summary>
    /// 评估上下文（可选）：用于在 NSGA-II 中设置“当前代数/个体编号”，
    /// 让 NSGAWiring 可以打印更精确的日志。
    ///
    /// 使用方法（可选）：
    ///   // 在 NSGAII.cs 的评估处（调用 Evaluate 之前）设置：
    ///   NSGAEvalContext.Set(genIndex: gen + 1, individualIndex: k + 1);
    ///   // 完成后可清理：
    ///   NSGAEvalContext.Clear();
    ///
    /// 若不设置，本文件会回退到“评估自增编号”的日志，不影响运行。
    /// </summary>
    public static class NSGAEvalContext
    {
        private static readonly AsyncLocal<int?> _generation = new();
        private static readonly AsyncLocal<int?> _individual = new();

        public static void Set(int? genIndex, int? individualIndex)
        {
            _generation.Value = genIndex;
            _individual.Value = individualIndex;
        }

        public static (int? gen, int? ind) Get() => (_generation.Value, _individual.Value);

        public static void Clear()
        {
            _generation.Value = null;
            _individual.Value = null;
        }
    }

    /// <summary>
    /// 将 GrowthSystem 与 LightingSimulator 串起来，提供 Evaluate(genes) 给 NSGA-II 使用。
    ///
    /// 基因布局：
    ///   genes[0..3]   → 4 个 repeller 因子，范围 [0.01, 5.0]
    ///   genes[4..403] → 400 个逐点位移（沿 -Y 法向），范围 [0, 100]
    ///
    /// 目标（统一最小化）：
    ///   f0 = 夏季光照小时（越小越好）
    ///   f1 = -冬季光照小时（冬季越多越好 → 取负）
    ///
    /// 新增：评估计时与进度输出
    ///   - 若未设置 NSGAEvalContext：打印 [评估 #K] 用时 XXX ms
    ///   - 若已设置 NSGAEvalContext：打印 [第 G 代 | 个体 I] 用时 XXX ms
    /// </summary>
    public static class NSGAWiring
    {
        // ======= GrowthSystem 默认参数 =======
        public static int    MaxPointCount = 200;
        public static int    MaxIterCount  = 200;
        public static double BaseDist      = 75.0;

        // ======= 几何参数 =======
        public static float  ExtrudeDepth  = 100f; // 沿 -Y 挤出的深度

        // ======= 光照模拟参数 =======
        public static DateOnly SummerDate  = new DateOnly(2025, 6, 21);
        public static DateOnly WinterDate  = new DateOnly(2025, 12, 21);
        public static TimeOnly StartTime   = new TimeOnly(8, 0);
        public static TimeOnly EndTime     = new TimeOnly(16, 0);
        public static TimeSpan Interval    = TimeSpan.FromHours(2);

        public static double RoomWidth     = 1000.0;
        public static double RoomDepth     = 1000.0;
        public static double GridSize      = 10.0;

        /// <summary>
        /// true：用平均日照小时作为目标；false：改为总日照（需你在 LightingSimulator 中实现 ComputeTotalLightHours）
        /// </summary>
        public static bool   UseAverageLightHours = false;

        // 全局评估计数（用于无上下文时打印“评估 #”）
        private static int _globalEvalCounter = 0;

        // =====================================================================
        // 工厂方法：生成 Evaluate 回调
        // =====================================================================

        /// <summary>
        /// 用内存点集创建 Evaluate，避免每次读盘。
        /// </summary>
        public static Func<double[], double[]> MakeEvaluator(
            List<Vector3> startingPoints,
            List<Vector3> repellerPoints)
        {
            if (startingPoints == null || startingPoints.Count == 0)
                throw new ArgumentException("startingPoints 不能为空。");
            if (repellerPoints == null)
                throw new ArgumentException("repellerPoints 不能为空。");

            return (genes) => EvaluateOnceWithLogging(genes, startingPoints, repellerPoints);
        }

        /// <summary>
        /// 从文件路径懒加载（首次 Evaluate 时加载并缓存）。
        /// </summary>
        public static Func<double[], double[]> MakeEvaluatorFromFiles(
            string startingCsvPath,
            string repellersCsvPath)
        {
            List<Vector3>? starting = null;
            List<Vector3>? repellers = null;

            return (genes) =>
            {
                starting  ??= IOHelper.LoadPointsFromFile(startingCsvPath);
                repellers ??= IOHelper.LoadPointsFromFile(repellersCsvPath);
                return EvaluateOnceWithLogging(genes, starting, repellers);
            };
        }

        // =====================================================================
        // 包裹一层：计时 + 控制台输出
        // =====================================================================

        private static double[] EvaluateOnceWithLogging(
            double[] genes,
            List<Vector3> startingPoints,
            List<Vector3> repellerPoints)
        {
            var sw = Stopwatch.StartNew();

            // 读取可选上下文（代数/个体编号），如果没有则回退到全局计数
            var (gen, ind) = NSGAEvalContext.Get();
            int evalId = Interlocked.Increment(ref _globalEvalCounter);

            // 真正的一次评估
            var result = EvaluateOnce(genes, startingPoints, repellerPoints);

            sw.Stop();

            if (gen.HasValue && ind.HasValue)
            {
                Console.WriteLine($"[第 {gen.Value} 代 | 个体 {ind.Value}] 用时 {sw.ElapsedMilliseconds} ms");
            }
            else
            {
                Console.WriteLine($"[评估 #{evalId}] 用时 {sw.ElapsedMilliseconds} ms");
            }

            return result;
        }

        // =====================================================================
        // 单次评估主流程：基因 → 平面生长 → 转竖直 → 逐点 -Y 偏移 → 挤出 -Y → 夏/冬光照 → 目标
        // =====================================================================

        private static double[] EvaluateOnce(
            double[] genes,
            List<Vector3> startingPoints,
            List<Vector3> repellerPoints)
        {
            if (genes == null || genes.Length < 404)
                throw new ArgumentException("基因长度不足（需要 404：4 个 repeller 因子 + 400 个逐点位移）。");

            // 1) 基因拆分
            const int repellerCount = 4;
            const int offsetCount   = 400;

            var repellerFactors = new List<double>(repellerCount);
            for (int i = 0; i < repellerCount; i++) repellerFactors.Add(genes[i]);

            var offsets = new double[offsetCount];
            Array.Copy(genes, repellerCount, offsets, 0, offsetCount);

            // 2) 平面生长（每次评估都重新计算）
            var growth = new GrowthSystem();
            var flatCurve = growth.Run(
                starting:        startingPoints,
                repellers:       repellerPoints,
                repellerFactors: repellerFactors,
                maxPointCount:   MaxPointCount,
                maxIterCount:    MaxIterCount,
                baseDist:        BaseDist
            );

            // 3) 转垂直（与你 Program.cs 一致）：(x, y, 0) → (x, 0, z=y)
            var verticalCrv = ToVerticalXZ(flatCurve);

            // 4) 逐点沿 -Y 偏移（只改前 N 个点）
            ApplyOffsetsMinusY(verticalCrv, offsets);

            // 5) 得到 extruded
            var extrudedCrv = verticalCrv;

            // 6) 夏 / 冬 光照模拟
            double summerMetric = SimulateAndGetMetric(verticalCrv, extrudedCrv, SummerDate);
            double winterMetric = SimulateAndGetMetric(verticalCrv, extrudedCrv, WinterDate);

            // 7) 统一最小化方向
            return new[] { summerMetric, -winterMetric };
        }

        // =====================================================================
        // 几何辅助
        // =====================================================================

        /// <summary>(x, y, 0) → (x, 0, z=y)</summary>
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

        private static double GetLightMetric(LightingSimulator sim)
        {
            if (UseAverageLightHours)
            {
                // 需要在 LightingSimulator 中提供：
                // public double ComputeAverageLightHours()
                return sim.GetAverageLightHours();
            }
            else
            {
                // 如果你实现了总计：
                // return sim.ComputeTotalLightHours();
                // 默认仍返回平均，避免编译错误：
                return sim.GetTotalLightHours();
            }
        }
    }
}
