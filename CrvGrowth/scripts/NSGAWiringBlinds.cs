// File: CrvGrowth/NSGAWiringBlinds.cs
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Numerics;
using System.Threading;

namespace CrvGrowth
{
    /// <summary>
    /// 专用于“百叶窗”几何的 NSGA-II 接线程序。
    ///
    /// 基因布局（总计 201）：
    ///   genes[0]        → 叶片个数 bladeCount（实值；实际使用时 round 后 clamp 到 [2, MaxBladeCount]）
    ///   genes[1..100]   → 最多 100 个叶片的原始拉伸长度候选值（建议范围 [0, UnitSize]）
    ///   genes[101..200] → 最多 100 个叶片的原始旋转角候选值（建议范围 [-85, 85]）
    ///
    /// 实际使用：
    ///   - bladeCount 确定后，只使用前 bladeCount 个长度 / 角度；
    ///   - 每个长度会被 clamp 到 [0, UnitSize / (bladeCount - 1)]；
    ///   - 每个角度会被 clamp 到 [-85, 85]。
    ///
    /// 目标（统一最小化）：
    ///   f0 = 夏季光照小时（越小越好）
    ///   f1 = -冬季光照小时（冬季越多越好 → 取负）
    /// </summary>
    public static class NSGAWiringBlinds
    {
        // ======= 几何参数 =======
        /// <summary>百叶单元的宽度 / 高度（mm），应与 BlindsGenerator 的 unitSize 一致。</summary>
        public static float UnitSize = 1000f;

        // ======= 光照模拟参数 =======
        public static DateOnly SummerDate  = new DateOnly(2025, 6, 21);
        public static DateOnly WinterDate  = new DateOnly(2025, 12, 21);
        public static TimeOnly StartTime   = new TimeOnly(9, 0);
        public static TimeOnly EndTime     = new TimeOnly(15, 0);
        public static TimeSpan Interval    = TimeSpan.FromHours(1);

        public static double RoomWidth     = 4000.0;
        public static double RoomDepth     = 3000.0;
        public static double GridSize      = 10.0;

        /// <summary>true：用平均日照小时作为目标；false：总日照</summary>
        public static bool   UseAverageLightHours = false;

        // 全局评估计数（用于无上下文时打印“评估 #”）
        private static int _globalEvalCounter = 0;

        // ======= 基因池布局常量 =======
        public const int MaxBladeCount = 100;

        private const int _BladeCountIndex   = 0;
        private const int _LengthStartIndex  = 1;
        private const int _LengthGeneCount   = MaxBladeCount;      // 1..100
        private const int _AngleStartIndex   = _LengthStartIndex + _LengthGeneCount; // 101
        private const int _AngleGeneCount    = MaxBladeCount;      // 101..200
        private const int _GeneLen           = 1 + _LengthGeneCount + _AngleGeneCount; // 201

        // =====================================================================
        // 工厂方法：生成 Evaluate 回调
        // =====================================================================

        /// <summary>
        /// 使用 NOAA / SPA 等内部路径在 Evaluate 内部计算太阳向量。
        /// </summary>
        public static Func<double[], double[]> MakeEvaluator()
        {
            return (genes) => EvaluateOnceWithLogging(genes);
        }

        /// <summary>
        /// 使用预计算好的太阳向量路径，避免 Evaluate 内重复天文计算。
        /// summerToSuns / winterToSuns 的含义与原 NSGAWiring 中一致。
        /// </summary>
        public static Func<double[], double[]> MakeEvaluatorUsingVectors(
            Vector3[] summerToSuns,
            Vector3[] winterToSuns)
        {
            if (summerToSuns == null || winterToSuns == null)
                throw new ArgumentException("太阳向量数组不能为空。");

            return (genes) => EvaluateOnceWithLogging_UsingVectors(genes, summerToSuns, winterToSuns);
        }

        // =====================================================================
        // 包裹一层：计时 + 控制台输出（NOAA / SPA 路径）
        // =====================================================================

        private static double[] EvaluateOnceWithLogging(double[] genes)
        {
            var sw = Stopwatch.StartNew();

            var (gen, ind) = NSGAEvalContext.Get(); // 复用原来的上下文
            int evalId = Interlocked.Increment(ref _globalEvalCounter);

            var result = EvaluateOnce(genes);

            sw.Stop();

            if (gen.HasValue && ind.HasValue)
                Console.WriteLine($"[百叶 | 第 {gen.Value} 代 | 个体 {ind.Value}] 用时 {sw.ElapsedMilliseconds} ms");
            else
                Console.WriteLine($"[百叶 | 评估 #{evalId}] 用时 {sw.ElapsedMilliseconds} ms");

            return result;
        }

        // =====================================================================
        // 包裹一层：计时 + 控制台输出（预计算太阳向量路径）
        // =====================================================================

        private static double[] EvaluateOnceWithLogging_UsingVectors(
            double[] genes,
            Vector3[] summerToSuns,
            Vector3[] winterToSuns)
        {
            var sw = Stopwatch.StartNew();

            var (gen, ind) = NSGAEvalContext.Get();
            int evalId = Interlocked.Increment(ref _globalEvalCounter);

            var result = EvaluateOnce_UsingVectors(genes, summerToSuns, winterToSuns);

            sw.Stop();

            if (gen.HasValue && ind.HasValue)
                Console.WriteLine($"[百叶 | 第 {gen.Value} 代 | 个体 {ind.Value}] 用时 {sw.ElapsedMilliseconds} ms");
            else
                Console.WriteLine($"[百叶 | 评估 #{evalId}] 用时 {sw.ElapsedMilliseconds} ms");

            return result;
        }

        // =====================================================================
        // 单次评估主流程（NOAA / SPA 路径）
        // =====================================================================

        private static double[] EvaluateOnce(double[] genes)
        {
            if (genes == null || genes.Length < _GeneLen)
                throw new ArgumentException($"基因长度不足（需要至少 {_GeneLen}：1 + 100 + 100）。");

            // 1) 解码 bladeCount
            int bladeCount = DecodeBladeCount(genes[_BladeCountIndex]);

            // 2) 使用基因生成百叶几何
            BlindsGenerator.GenerateFromGenes(
                bladeCount:       bladeCount,
                unitSize:         UnitSize,
                genes:            genes,
                lengthStartIndex: _LengthStartIndex,
                angleStartIndex:  _AngleStartIndex,
                out var verticalCrv,
                out var extrudedCrv
            );

            // 3) 夏 / 冬 光照模拟（百叶窗模式）
            double summerMetric = SimulateAndGetMetric(verticalCrv, extrudedCrv, SummerDate);
            double winterMetric = SimulateAndGetMetric(verticalCrv, extrudedCrv, WinterDate);

            // 4) 统一最小化方向
            return new[] { summerMetric, -winterMetric };
        }

        // =====================================================================
        // 单次评估主流程（预计算太阳向量路径）
        // =====================================================================

        private static double[] EvaluateOnce_UsingVectors(
            double[] genes,
            Vector3[] summerToSuns,
            Vector3[] winterToSuns)
        {
            if (genes == null || genes.Length < _GeneLen)
                throw new ArgumentException($"基因长度不足（需要至少 {_GeneLen}：1 + 100 + 100）。");

            int bladeCount = DecodeBladeCount(genes[_BladeCountIndex]);

            BlindsGenerator.GenerateFromGenes(
                bladeCount:       bladeCount,
                unitSize:         UnitSize,
                genes:            genes,
                lengthStartIndex: _LengthStartIndex,
                angleStartIndex:  _AngleStartIndex,
                out var verticalCrv,
                out var extrudedCrv
            );

            double summerMetric = SimulateAndGetMetric_WithVectors(verticalCrv, extrudedCrv, summerToSuns);
            double winterMetric = SimulateAndGetMetric_WithVectors(verticalCrv, extrudedCrv, winterToSuns);

            return new[] { summerMetric, -winterMetric };
        }

        // =====================================================================
        // 解码辅助
        // =====================================================================

        private static int DecodeBladeCount(double raw)
        {
            int bc = (int)Math.Round(raw);
            if (bc < 2) bc = 2;
            if (bc > MaxBladeCount) bc = MaxBladeCount;
            return bc;
        }

        // =====================================================================
        // 光照辅助（NOAA / SPA 路径 → 使用 RunSimulationBlinds）
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
                gridSize:      GridSize,
                isClosed:      false,  // 百叶为一片片，不闭合
                enablePeriodicTiling: true,
                mirrorOffset:  UnitSize // 一般与单元宽度一致
            );

            sim.RunSimulationBlinds();
            return GetLightMetric(sim);
        }

        // =====================================================================
        // 光照辅助（预计算太阳向量路径 → 使用 RunWithSunVectorsBlinds）
        // =====================================================================

        private static double SimulateAndGetMetric_WithVectors(
            List<Vector3> verticalCrv,
            List<Vector3> extrudedCrv,
            Vector3[] toSuns)
        {
            var sim = new LightingSimulator(
                verticalCurve: verticalCrv,
                extrudedCurve: extrudedCrv,
                date:          SummerDate,  // 占位，不在 vector 模式中实际使用
                startTime:     StartTime,
                endTime:       EndTime,
                interval:      Interval,
                roomWidth:     RoomWidth,
                roomDepth:     RoomDepth,
                gridSize:      GridSize,
                isClosed:      false,
                enablePeriodicTiling: true,
                mirrorOffset:  UnitSize
            );

            sim.RunWithSunVectorsBlinds(toSuns);
            return GetLightMetric(sim);
        }

        private static double GetLightMetric(LightingSimulator sim)
        {
            return UseAverageLightHours ? sim.GetAverageLightHours() : sim.GetTotalLightHours();
        }
    }
}
