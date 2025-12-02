// File: CrvGrowth/BlindsGenerator.cs
using System;
using System.Collections.Generic;
using System.Numerics;

namespace CrvGrowth
{
    public static class BlindsGenerator
    {
        /// <summary>
        /// 生成百叶窗的基础点集（未拉伸、未旋转）。
        /// 坐标位于 xz 平面，y = 0。
        /// 顺序：L0, R0, L1, R1, ..., L(n-1), R(n-1)
        /// </summary>
        public static List<Vector3> GenerateBasePoints(int bladeCount, float unitSize = 1000f)
        {
            if (bladeCount <= 1)
                throw new ArgumentException("叶片个数必须大于 1。", nameof(bladeCount));

            var points = new List<Vector3>(bladeCount * 2);

            float segmentCount = bladeCount - 1;
            float dz = unitSize / segmentCount;

            for (int i = 0; i < bladeCount; i++)
            {
                float z = dz * i;

                // 左边界点 (x = 0)
                points.Add(new Vector3(0f, 0f, z));

                // 右边界点 (x = unitSize)
                points.Add(new Vector3(unitSize, 0f, z));
            }

            return points;
        }

        /// <summary>
        /// 绕任意轴旋转一个点。
        /// axisPoint: 轴上一点
        /// axisDirNormalized: 归一化的轴方向向量
        /// angleRad: 旋转角度（弧度）
        /// </summary>
        private static Vector3 RotateAroundAxis(
            Vector3 point,
            Vector3 axisPoint,
            Vector3 axisDirNormalized,
            float angleRad)
        {
            if (axisDirNormalized.LengthSquared() < 1e-8f)
                return point;

            var v = axisDirNormalized;
            var p = point - axisPoint;

            float cos = MathF.Cos(angleRad);
            float sin = MathF.Sin(angleRad);

            var pParallel = Vector3.Dot(p, v) * v;
            var pPerp = p - pParallel;
            var pPerpRot = pPerp * cos + Vector3.Cross(v, pPerp) * sin;

            var rotated = axisPoint + pParallel + pPerpRot;
            return rotated;
        }

        /// <summary>
        /// 示例：为每片叶片随机生成拉伸长度和角度，
        /// length[i] ∈ [0, unitSize / (bladeCount - 1)]
        /// angle[i]  ∈ [-85°, +85°]
        ///
        /// 输出：
        ///  basePoints     : 基础铰接线点（2 * bladeCount）
        ///  extrudedPoints : 沿 -Y 拉伸后的点（与 basePoints 同长度、同顺序）
        ///  rotatedPoints  : 在拉伸基础上绕各自轴旋转后的点（与 basePoints 同长度、同顺序）
        /// </summary>
        public static void GenerateWithExtrudeAndRotateRandomExample(
            int bladeCount,
            float unitSize,
            out List<Vector3> basePoints,
            out List<Vector3> extrudedPoints,
            out List<Vector3> rotatedPoints,
            out float[] lengths,
            out float[] anglesDeg)
        {
            if (bladeCount <= 1)
                throw new ArgumentException("叶片个数必须大于 1。", nameof(bladeCount));

            basePoints = GenerateBasePoints(bladeCount, unitSize);

            int pointCount = basePoints.Count;             // 应为 2 * bladeCount
            int expectedPointCount = bladeCount * 2;

            if (pointCount != expectedPointCount)
                throw new InvalidOperationException(
                    $"基础点数量异常: {pointCount} != {expectedPointCount}");

            extrudedPoints = new List<Vector3>(pointCount);
            rotatedPoints  = new List<Vector3>(pointCount);

            lengths   = new float[bladeCount];
            anglesDeg = new float[bladeCount];

            var rand = new Random();

            float maxLength = unitSize / (bladeCount - 1);
            const float EPS = 1e-8f;

            for (int i = 0; i < bladeCount; i++)
            {
                // 随机生成当前叶片的拉伸长度和角度
                float length = (float)(rand.NextDouble() * maxLength);              // [0, maxLength]
                float angle  = (float)(-85.0 + rand.NextDouble() * 170.0);         // [-85, +85]

                lengths[i]   = length;
                anglesDeg[i] = angle;

                int idxLeft  = 2 * i;
                int idxRight = 2 * i + 1;

                var baseLeft  = basePoints[idxLeft];
                var baseRight = basePoints[idxRight];

                // 1) 沿 -Y 方向拉伸
                var offset = new Vector3(0f, -length, 0f);
                var extrudedLeft  = baseLeft  + offset;
                var extrudedRight = baseRight + offset;

                extrudedPoints.Add(extrudedLeft);
                extrudedPoints.Add(extrudedRight);

                // 2) 围绕“原始左右点连线”为轴进行旋转
                var axisDir = baseRight - baseLeft;
                Vector3 rotatedLeft;
                Vector3 rotatedRight;

                if (axisDir.LengthSquared() > EPS)
                {
                    axisDir = Vector3.Normalize(axisDir);
                    float angleRad = angle * (MathF.PI / 180f);

                    rotatedLeft  = RotateAroundAxis(extrudedLeft,  baseLeft, axisDir, angleRad);
                    rotatedRight = RotateAroundAxis(extrudedRight, baseLeft, axisDir, angleRad);
                }
                else
                {
                    // 极端情况下左右点重合，则不旋转
                    rotatedLeft  = extrudedLeft;
                    rotatedRight = extrudedRight;
                }

                // 无论走哪个分支，这里统一加入两点，避免漏点
                rotatedPoints.Add(rotatedLeft);
                rotatedPoints.Add(rotatedRight);
            }

            // 数量一致性检查
            if (extrudedPoints.Count != expectedPointCount)
            {
                throw new InvalidOperationException(
                    $"Extruded 点数量异常: {extrudedPoints.Count} != {expectedPointCount}");
            }

            if (rotatedPoints.Count != expectedPointCount)
            {
                throw new InvalidOperationException(
                    $"Rotated 点数量异常: {rotatedPoints.Count} != {expectedPointCount}");
            }
        }

        /// <summary>
        /// 从 NSGA 基因向量中生成百叶几何：
        /// - bladeCount        : 当前使用的叶片个数（外部已 clamp 到 [2, max]）
        /// - unitSize          : 窗口单元宽度 / 高度（mm）
        /// - genes             : NSGA 基因向量
        /// - lengthStartIndex  : genes 中长度基因起始索引（例如 1）
        /// - angleStartIndex   : genes 中角度基因起始索引（例如 101）
        ///
        /// 输出：
        ///   basePoints   → 作为 verticalCurve（位于 xz 平面，y = 0）
        ///   extrudedEdge → 作为 extrudedCurve（沿 -Y 拉伸并绕轴旋转后的外缘）
        ///
        /// 实际使用：
        ///   - 每片叶片的最大拉伸长度 = unitSize / (bladeCount - 1)，
        ///   - 真实长度 clamp 到 [0, maxLength]，
        ///   - 角度 clamp 到 [-85, 85]。
        /// </summary>
        public static void GenerateFromGenes(
            int bladeCount,
            float unitSize,
            double[] genes,
            int lengthStartIndex,
            int angleStartIndex,
            out List<Vector3> basePoints,
            out List<Vector3> extrudedEdge)
        {
            if (bladeCount <= 1)
                throw new ArgumentException("bladeCount 必须大于 1。", nameof(bladeCount));
            if (genes == null)
                throw new ArgumentNullException(nameof(genes));

            basePoints = GenerateBasePoints(bladeCount, unitSize);

            int pointCount = basePoints.Count; // = 2 * bladeCount
            int expectedPointCount = bladeCount * 2;
            if (pointCount != expectedPointCount)
                throw new InvalidOperationException(
                    $"基础点数量异常: {pointCount} != {expectedPointCount}");

            extrudedEdge = new List<Vector3>(pointCount);

            double maxLength = unitSize / (bladeCount - 1);
            const float EPS = 1e-8f;

            for (int i = 0; i < bladeCount; i++)
            {
                int idxLeft  = 2 * i;
                int idxRight = 2 * i + 1;

                var baseLeft  = basePoints[idxLeft];
                var baseRight = basePoints[idxRight];

                int lenIndex = lengthStartIndex + i;
                int angIndex = angleStartIndex + i;

                if (lenIndex >= genes.Length || angIndex >= genes.Length)
                    throw new ArgumentException("基因向量长度不足以支撑 bladeCount 对应的长度/角度。");

                double rawLength = genes[lenIndex];
                double rawAngle  = genes[angIndex];

                float length   = (float)Math.Clamp(rawLength, 0.0, maxLength);
                float angleDeg = (float)Math.Clamp(rawAngle, -85.0, 85.0);

                var offset = new Vector3(0f, -length, 0f);
                var extrudedLeft  = baseLeft  + offset;
                var extrudedRight = baseRight + offset;

                var axisDir = baseRight - baseLeft;
                Vector3 rotatedLeft;
                Vector3 rotatedRight;

                if (axisDir.LengthSquared() > EPS)
                {
                    axisDir = Vector3.Normalize(axisDir);
                    float angleRad = angleDeg * (MathF.PI / 180f);

                    rotatedLeft  = RotateAroundAxis(extrudedLeft,  baseLeft, axisDir, angleRad);
                    rotatedRight = RotateAroundAxis(extrudedRight, baseLeft, axisDir, angleRad);
                }
                else
                {
                    rotatedLeft  = extrudedLeft;
                    rotatedRight = extrudedRight;
                }

                extrudedEdge.Add(rotatedLeft);
                extrudedEdge.Add(rotatedRight);
            }

            if (extrudedEdge.Count != expectedPointCount)
            {
                throw new InvalidOperationException(
                    $"extrudedEdge 点数量异常: {extrudedEdge.Count} != {expectedPointCount}");
            }
        }
    }
}
