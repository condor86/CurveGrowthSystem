// File: CrvGrowth/BlindsGenerator.cs
using System;
using System.Collections.Generic;
using System.Numerics;

namespace CrvGrowth
{
    /// <summary>
    /// 百叶窗几何生成：
    /// - 基础点：位于 XZ 平面，左边界 x=0，右边界 x=unitSize，Z 从 0 到 unitSize 等分。
    /// - 点序：L0, R0, L1, R1, ..., L_{n-1}, R_{n-1}
    /// - 从基因解码：先沿 -Y 拉伸，再绕各自“左右边界连线”旋转。
    /// </summary>
    public static class BlindsGenerator
    {
        public const float UnitSizeDefault = 1000f;

        /// <summary>
        /// 生成基础百叶点集（仅按照 bladeCount 等分 Z，尚未拉伸/旋转）。
        /// 左下角在原点 (0,0,0)，右上角在 (unitSize, 0, unitSize)。
        /// 点序：L0, R0, L1, R1, ..., L_{bladeCount-1}, R_{bladeCount-1}
        /// </summary>
        public static List<Vector3> GenerateBaseBlindsPoints(int bladeCount, float unitSize)
        {
            if (bladeCount < 2)
                throw new ArgumentOutOfRangeException(nameof(bladeCount), "bladeCount 必须 >= 2。");
            if (unitSize <= 0)
                throw new ArgumentOutOfRangeException(nameof(unitSize), "unitSize 必须 > 0。");

            var points = new List<Vector3>(bladeCount * 2);
            float step = unitSize / (bladeCount - 1);

            for (int i = 0; i < bladeCount; i++)
            {
                float z = step * i;
                // 左点 (x=0, y=0, z)
                points.Add(new Vector3(0f, 0f, z));
                // 右点 (x=unitSize, y=0, z)
                points.Add(new Vector3(unitSize, 0f, z));
            }

            return points;
        }

        /// <summary>
        /// 用于单独测试百叶生成逻辑的随机版本。
        /// 每片叶片：
        /// - 长度 ∈ [0, unitSize/(bladeCount-1)]
        /// - 角度 ∈ [-85°, 85°]
        /// </summary>
        public static void GenerateRandomBlinds(
            int bladeCount,
            float unitSize,
            out List<Vector3> basePoints,
            out List<Vector3> extrudedEdge,
            Random? rng = null)
        {
            rng ??= new Random();

            basePoints = GenerateBaseBlindsPoints(bladeCount, unitSize);
            extrudedEdge = new List<Vector3>(basePoints.Count);

            double maxLength = unitSize / (bladeCount - 1);
            const float EPS = 1e-8f;

            for (int i = 0; i < bladeCount; i++)
            {
                int idxLeft  = 2 * i;
                int idxRight = 2 * i + 1;

                var baseLeft  = basePoints[idxLeft];
                var baseRight = basePoints[idxRight];

                float length   = (float)(rng.NextDouble() * maxLength);
                float angleDeg = (float)(rng.NextDouble() * 170.0 - 85.0); // [-85, +85]

                var offset = new Vector3(0f, -length, 0f);
                var extrudedLeft  = baseLeft  + offset;
                var extrudedRight = baseRight + offset;

                var axisDir = baseRight - baseLeft;
                if (axisDir.LengthSquared() < EPS)
                {
                    // 退化：不旋转，直接使用拉伸结果
                    extrudedEdge.Add(extrudedLeft);
                    extrudedEdge.Add(extrudedRight);
                    continue;
                }

                var axisUnit = Vector3.Normalize(axisDir);
                float angleRad = angleDeg * (float)(Math.PI / 180.0);

                var rotatedLeft  = RotateAroundAxisThroughPoint(extrudedLeft,  baseLeft, axisUnit, angleRad);
                var rotatedRight = RotateAroundAxisThroughPoint(extrudedRight, baseLeft, axisUnit, angleRad);

                extrudedEdge.Add(rotatedLeft);
                extrudedEdge.Add(rotatedRight);
            }
        }

        /// <summary>
        /// 1) 先根据 bladeCount 计算 maxLength = unitSize / (bladeCount - 1)；
        /// 2) genes[lenIndex] 视为在 [0, unitSize] 上的随机数；
        /// 3) 在此处线性映射到 [0, maxLength]；
        /// 4) 只取前 bladeCount 个长度和角度生成几何。
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
            if (genes == null)
                throw new ArgumentNullException(nameof(genes));
            if (bladeCount < 2)
                throw new ArgumentOutOfRangeException(nameof(bladeCount), "bladeCount 必须 >= 2。");
            if (unitSize <= 0)
                throw new ArgumentOutOfRangeException(nameof(unitSize), "unitSize 必须 > 0。");

            basePoints   = GenerateBaseBlindsPoints(bladeCount, unitSize);
            extrudedEdge = new List<Vector3>(basePoints.Count);

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
                    throw new ArgumentException("基因长度不足以支撑 bladeCount 对应的长度/角度。");

                double rawLength = genes[lenIndex];
                double rawAngle  = genes[angIndex];
                
                // rawLength 先视为 [0, unitSize] 上的随机值，
                // 再缩放到 [0, maxLength]。
                double t = rawLength / unitSize;   // 变成 0..1 比例
                t = Math.Clamp(t, 0.0, 1.0);       // 防越界
                float length   = (float)(t * maxLength);
                float angleDeg = (float)Math.Clamp(rawAngle, -85.0, 85.0);

                // 先沿 -Y 拉伸
                var offset = new Vector3(0f, -length, 0f);
                var extrudedLeft  = baseLeft  + offset;
                var extrudedRight = baseRight + offset;

                // 再绕“原始左右两点连线”为轴旋转
                var axisDir = baseRight - baseLeft;
                if (axisDir.LengthSquared() < EPS)
                {
                    // 极端退化：不旋转，只保留拉伸后的结果，依然保留两点
                    extrudedEdge.Add(extrudedLeft);
                    extrudedEdge.Add(extrudedRight);
                    continue;
                }

                var axisUnit = Vector3.Normalize(axisDir);
                float angleRad = angleDeg * (float)(Math.PI / 180.0);

                var rotatedLeft  = RotateAroundAxisThroughPoint(extrudedLeft,  baseLeft, axisUnit, angleRad);
                var rotatedRight = RotateAroundAxisThroughPoint(extrudedRight, baseLeft, axisUnit, angleRad);

                extrudedEdge.Add(rotatedLeft);
                extrudedEdge.Add(rotatedRight);
            }
        }

        /// <summary>
        /// 先把点移到以 axisOrigin 为原点的坐标系，再绕单位轴 axisUnit 旋转 angleRad，然后移回去。
        /// </summary>
        private static Vector3 RotateAroundAxisThroughPoint(
            in Vector3 point,
            in Vector3 axisOrigin,
            in Vector3 axisUnit,
            float angleRad)
        {
            var v    = point - axisOrigin;
            var vRot = RotateAroundAxis(v, axisUnit, angleRad);
            return axisOrigin + vRot;
        }

        /// <summary>Rodrigues 公式：绕单位轴 axisUnit 旋转 angleRad。</summary>
        private static Vector3 RotateAroundAxis(
            in Vector3 v,
            in Vector3 axisUnit,
            float angleRad)
        {
            float c = MathF.Cos(angleRad);
            float s = MathF.Sin(angleRad);

            return v * c
                 + Vector3.Cross(axisUnit, v) * s
                 + axisUnit * Vector3.Dot(axisUnit, v) * (1 - c);
        }
    }
}
