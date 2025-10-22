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
    /// 评估上下文（可选）：用于在 NSGA-II 中设置“当前代数/个体编号”，便于日志。
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
    /// 新基因布局（总计 804）：
    ///   genes[0..3]     → 4 个 repeller 因子，范围 [0.01, 5.0]
    ///   genes[4..403]   → 400 个逐点位移（沿 -Y 法向），范围 [50, 100]（与 Program 保持一致）
    ///   genes[404..803] → 400 个逐点“局部平面旋转”角度（单位：度），范围 [-25, 25]
    /// 
    /// 目标（统一最小化）：
    ///   f0 = 夏季光照小时（越小越好）
    ///   f1 = -冬季光照小时（冬季越多越好 → 取负）
    /// 
    /// 新增：可选“预计算太阳向量”通道，避免在热循环内重复天文计算
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
        public static TimeOnly StartTime   = new TimeOnly(9, 0);
        public static TimeOnly EndTime     = new TimeOnly(15, 0);
        public static TimeSpan Interval    = TimeSpan.FromHours(1);

        public static double RoomWidth     = 1000.0;
        public static double RoomDepth     = 1000.0;
        public static double GridSize      = 10.0;

        /// <summary>
        /// true：用平均日照小时作为目标；false：改为总日照
        /// </summary>
        public static bool   UseAverageLightHours = false;

        // 全局评估计数（用于无上下文时打印“评估 #”）
        private static int _globalEvalCounter = 0;

        // 索引常量
        private const int _RepellerCount = 4;
        private const int _OffsetCount   = 400;
        private const int _AngleCount    = 400;
        private const int _GeneLen       = _RepellerCount + _OffsetCount + _AngleCount; // 804

        // =====================================================================
        // 工厂方法：生成 Evaluate 回调
        // =====================================================================

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

        /// <summary>
        /// ✅ 预计算太阳向量路径（推荐）
        /// </summary>
        public static Func<double[], double[]> MakeEvaluator(
            List<Vector3> startingPoints,
            List<Vector3> repellerPoints,
            Vector3[] summerToSuns,
            Vector3[] winterToSuns)
        {
            if (startingPoints == null || startingPoints.Count == 0)
                throw new ArgumentException("startingPoints 不能为空。");
            if (repellerPoints == null)
                throw new ArgumentException("repellerPoints 不能为空。");
            if (summerToSuns == null || winterToSuns == null)
                throw new ArgumentException("太阳向量数组不能为空。");

            return (genes) => EvaluateOnceWithLogging_UsingVectors(
                genes, startingPoints, repellerPoints, summerToSuns, winterToSuns);
        }

        // =====================================================================
        // 包裹一层：计时 + 控制台输出（原始 NOAA 路径）
        // =====================================================================

        private static double[] EvaluateOnceWithLogging(
            double[] genes,
            List<Vector3> startingPoints,
            List<Vector3> repellerPoints)
        {
            var sw = Stopwatch.StartNew();

            var (gen, ind) = NSGAEvalContext.Get();
            int evalId = Interlocked.Increment(ref _globalEvalCounter);

            var result = EvaluateOnce(genes, startingPoints, repellerPoints);

            sw.Stop();

            if (gen.HasValue && ind.HasValue)
                Console.WriteLine($"[第 {gen.Value} 代 | 个体 {ind.Value}] 用时 {sw.ElapsedMilliseconds} ms");
            else
                Console.WriteLine($"[评估 #{evalId}] 用时 {sw.ElapsedMilliseconds} ms");

            return result;
        }

        // =====================================================================
        // 包裹一层：计时 + 控制台输出（预计算向量路径）
        // =====================================================================

        private static double[] EvaluateOnceWithLogging_UsingVectors(
            double[] genes,
            List<Vector3> startingPoints,
            List<Vector3> repellerPoints,
            Vector3[] summerToSuns,
            Vector3[] winterToSuns)
        {
            var sw = Stopwatch.StartNew();

            var (gen, ind) = NSGAEvalContext.Get();
            int evalId = Interlocked.Increment(ref _globalEvalCounter);

            var result = EvaluateOnce_UsingVectors(
                genes, startingPoints, repellerPoints, summerToSuns, winterToSuns);

            sw.Stop();

            if (gen.HasValue && ind.HasValue)
                Console.WriteLine($"[第 {gen.Value} 代 | 个体 {ind.Value}] 用时 {sw.ElapsedMilliseconds} ms");
            else
                Console.WriteLine($"[评估 #{evalId}] 用时 {sw.ElapsedMilliseconds} ms");

            return result;
        }

        // =====================================================================
        // 单次评估主流程（原始 NOAA 路径）
        // =====================================================================

        private static double[] EvaluateOnce(
            double[] genes,
            List<Vector3> startingPoints,
            List<Vector3> repellerPoints)
        {
            if (genes == null || genes.Length < _GeneLen)
                throw new ArgumentException($"基因长度不足（需要 {_GeneLen}：4 + 400 + 400）。");

            // 1) 基因拆分
            var repellerFactors = new List<double>(_RepellerCount);
            for (int i = 0; i < _RepellerCount; i++) repellerFactors.Add(genes[i]);

            var offsets = new double[_OffsetCount];
            Array.Copy(genes, _RepellerCount, offsets, 0, _OffsetCount);

            var anglesDeg = new double[_AngleCount];
            Array.Copy(genes, _RepellerCount + _OffsetCount, anglesDeg, 0, _AngleCount);

            // 2) 生形（平面）
            var growth = new GrowthSystem();
            var flatCurve = growth.Run(
                starting:        startingPoints,
                repellers:       repellerPoints,
                repellerFactors: repellerFactors,
                maxPointCount:   MaxPointCount,
                maxIterCount:    MaxIterCount,
                baseDist:        BaseDist
            );

            // 3) 转垂直
            var verticalCrv = ToVerticalXZ(flatCurve);

            // 4) 逐点 -Y 偏移（前 N）
            var extrudedCrv = new List<Vector3>(verticalCrv);
            int N = Math.Min(verticalCrv.Count, offsets.Length);
            for (int i = 0; i < N; i++)
            {
                var p = verticalCrv[i];
                extrudedCrv[i] = new Vector3(p.X, p.Y - (float)offsets[i], p.Z);
            }

            // 5) 新增：逐点局部平面旋转（以 Pn 为枢轴；平面由角平分线与 pnPn 线生成）
            ApplyLocalPlaneRotation(verticalCrv, extrudedCrv, anglesDeg);

            // 6) 夏 / 冬 光照模拟（NOAA 现算）
            double summerMetric = SimulateAndGetMetric(verticalCrv, extrudedCrv, SummerDate);
            double winterMetric = SimulateAndGetMetric(verticalCrv, extrudedCrv, WinterDate);

            // 7) 统一最小化方向
            return new[] { summerMetric, -winterMetric };
        }

        // =====================================================================
        // 单次评估主流程（预计算向量路径）
        // =====================================================================

        private static double[] EvaluateOnce_UsingVectors(
            double[] genes,
            List<Vector3> startingPoints,
            List<Vector3> repellerPoints,
            Vector3[] summerToSuns,
            Vector3[] winterToSuns)
        {
            if (genes == null || genes.Length < _GeneLen)
                throw new ArgumentException($"基因长度不足（需要 {_GeneLen}：4 + 400 + 400）。");

            var repellerFactors = new List<double>(_RepellerCount);
            for (int i = 0; i < _RepellerCount; i++) repellerFactors.Add(genes[i]);

            var offsets = new double[_OffsetCount];
            Array.Copy(genes, _RepellerCount, offsets, 0, _OffsetCount);

            var anglesDeg = new double[_AngleCount];
            Array.Copy(genes, _RepellerCount + _OffsetCount, anglesDeg, 0, _AngleCount);

            // 生形
            var growth = new GrowthSystem();
            var flatCurve = growth.Run(
                starting:        startingPoints,
                repellers:       repellerPoints,
                repellerFactors: repellerFactors,
                maxPointCount:   MaxPointCount,
                maxIterCount:    MaxIterCount,
                baseDist:        BaseDist
            );

            // 转垂直 + -Y 偏移
            var verticalCrv = ToVerticalXZ(flatCurve);
            var extrudedCrv = new List<Vector3>(verticalCrv);
            int N = Math.Min(verticalCrv.Count, offsets.Length);
            for (int i = 0; i < N; i++)
            {
                var p = verticalCrv[i];
                extrudedCrv[i] = new Vector3(p.X, p.Y - (float)offsets[i], p.Z);
            }

            // 新增：逐点局部平面旋转
            ApplyLocalPlaneRotation(verticalCrv, extrudedCrv, anglesDeg);

            // 光照：直接用预计算的向量数组
            double summerMetric = SimulateAndGetMetric_WithVectors(verticalCrv, extrudedCrv, summerToSuns);
            double winterMetric = SimulateAndGetMetric_WithVectors(verticalCrv, extrudedCrv, winterToSuns);

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

        /// <summary>
        /// 新增：逐点“局部平面旋转”
        /// - 枢轴：Pn = verticalCrv[i]
        /// - 初始向量：v = pn - Pn（pn 来自 extrudedCrv[i]）
        /// - 平面：由“角平分线 b̂（pn-1,pn,pn+1）”与“d̂ = normalize(Pn - pn)”张成
        /// - 旋转轴：n̂ = normalize( b̂ × d̂ )（平面法向）
        /// - 旋转角：anglesDeg[i]（度），右手法则；若与直觉相反，可在基因中取负或在此处统一反号
        /// </summary>
        public static void ApplyLocalPlaneRotation(
            List<Vector3> verticalCrv,
            List<Vector3> extrudedCrv,
            double[] anglesDeg)
        {
            if (verticalCrv == null || extrudedCrv == null || anglesDeg == null) return;

            int n = Math.Min(Math.Min(verticalCrv.Count, extrudedCrv.Count), anglesDeg.Length);
            if (n <= 0) return;

            const float EPS = 1e-8f;

            for (int i = 0; i < n; i++)
            {
                var Pn = verticalCrv[i];   // 枢轴
                var pn = extrudedCrv[i];   // 当前点（偏移后）
                var v  = pn - Pn;          // 旋转向量

                if (v.LengthSquared() < EPS) continue; // 与枢轴重合，无需旋转

                // d̂：从 pn 指向 Pn（与 v 反向）
                var d = Pn - pn;
                float dlen2 = d.LengthSquared();
                if (dlen2 < EPS) continue;
                var dHat = Vector3.Normalize(d);

                // 角平分线方向 b̂（基于 extrudedCrv 的局部几何）
                Vector3 bHat;
                bool hasL = (i - 1) >= 0;
                bool hasR = (i + 1) < extrudedCrv.Count;

                if (hasL && hasR)
                {
                    var d1 = extrudedCrv[i]   - extrudedCrv[i - 1];
                    var d2 = extrudedCrv[i + 1] - extrudedCrv[i];
                    if (d1.LengthSquared() > EPS) d1 = Vector3.Normalize(d1);
                    if (d2.LengthSquared() > EPS) d2 = Vector3.Normalize(d2);
                    var sum = d1 + d2;
                    if (sum.LengthSquared() < EPS)
                    {
                        // 180° 或极端情况：退化为单侧切向
                        bHat = (d2.LengthSquared() > EPS) ? d2 :
                               (d1.LengthSquared() > EPS) ? d1 : new Vector3(1, 0, 0);
                    }
                    else bHat = Vector3.Normalize(sum);
                }
                else if (hasR || hasL)
                {
                    var tan = hasR ? (extrudedCrv[i + 1] - extrudedCrv[i]) : (extrudedCrv[i] - extrudedCrv[i - 1]);
                    if (tan.LengthSquared() < EPS) continue;
                    bHat = Vector3.Normalize(tan);
                }
                else
                {
                    // 仅 1 点：随便取一条不与 d̂ 共线的方向
                    bHat = Math.Abs(dHat.X) < 0.9f ? Vector3.UnitX : Vector3.UnitY;
                }

                // 法向 n̂ = b̂ × d̂
                var axis = Vector3.Cross(bHat, dHat);
                if (axis.LengthSquared() < EPS)
                {
                    // 共线退化：取与 d̂ 不共线的全局轴做叉积
                    var fallback = Math.Abs(dHat.X) < 0.9f ? Vector3.UnitX : Vector3.UnitY;
                    axis = Vector3.Cross(fallback, dHat);
                    if (axis.LengthSquared() < EPS)
                        axis = Vector3.Cross(Vector3.UnitZ, dHat);
                    if (axis.LengthSquared() < EPS) continue; // 仍退化则跳过
                }
                axis = Vector3.Normalize(axis);

                // 旋转角（度→弧度）
                float theta = (float)(anglesDeg[i] * Math.PI / 180.0);

                // Rodrigues 旋转（绕 axis，通过 Pn）
                var vRot = RotateAroundAxis(v, axis, theta);

                extrudedCrv[i] = Pn + vRot;
            }
        }

        /// <summary>Rodrigues 旋转公式：绕单位轴 axis 旋转 angleRad</summary>
        private static Vector3 RotateAroundAxis(in Vector3 v, in Vector3 axis, float angleRad)
        {
            float c = MathF.Cos(angleRad);
            float s = MathF.Sin(angleRad);
            return v * c + Vector3.Cross(axis, v) * s + axis * Vector3.Dot(axis, v) * (1 - c);
        }

        // =====================================================================
        // 光照辅助（原始 NOAA 路径）
        // =====================================================================

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

        // =====================================================================
        // 光照辅助（预计算向量路径）
        // =====================================================================

        private static double SimulateAndGetMetric_WithVectors(
            List<Vector3> verticalCrv,
            List<Vector3> extrudedCrv,
            Vector3[] toSuns)
        {
            var sim = new LightingSimulator(
                verticalCurve: verticalCrv,
                extrudedCurve: extrudedCrv,
                date:          SummerDate,  // 占位，不在该路径中使用
                startTime:     StartTime,
                endTime:       EndTime,
                interval:      Interval,
                roomWidth:     RoomWidth,
                roomDepth:     RoomDepth,
                gridSize:      GridSize
            );
            sim.RunWithSunVectors(toSuns);
            return GetLightMetric(sim);
        }

        private static double GetLightMetric(LightingSimulator sim)
        {
            return UseAverageLightHours ? sim.GetAverageLightHours() : sim.GetTotalLightHours();
        }
    }
}
