using Supercluster.KDTree;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace CrvGrowth
{
    internal static class TopologyHelpers
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static int PrevIndex(int i, int n, bool closed)
            => closed ? (i - 1 + n) % n : Math.Max(i - 1, 0);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static int NextIndex(int i, int n, bool closed)
            => closed ? (i + 1) % n : Math.Min(i + 1, n - 1);

        /// <summary>
        /// 段枚举：当 closed==true 时返回 n 段（含 last→0），否则返回 n-1 段（0→1, ..., n-2→n-1）。
        /// </summary>
        public static IEnumerable<(int a, int b)> EnumerateSegments(int n, bool closed)
        {
            if (n < 2) yield break;

            if (closed)
            {
                for (int i = 0; i < n; i++)
                    yield return (i, (i + 1) % n);
            }
            else
            {
                for (int i = 0; i < n - 1; i++)
                    yield return (i, i + 1);
            }
        }
    }

    public class GrowthSystem
    {
        private readonly double _tileWidth = 1000.0;
        private readonly double _tileHeight = 1000.0;
        private readonly double _maxFactor = 1.5;
        private readonly double _maxEffectDist = 300.0;

        public List<Vector3> Run(
            List<Vector3> starting,
            List<Vector3> repellers,
            List<double> repellerFactors,
            int maxPointCount = 200,
            int maxIterCount = 200,
            double baseDist = 75.0,
            bool isClosed = true) // 新增：闭合拓扑开关（默认开启）
        {
            var centers = new List<Vector3>(starting);

            // 若希望闭合且起始点数过少，建议上游保障 >= 3；此处不强制改动以保持兼容
            for (int iter = 0; iter < maxIterCount; iter++)
            {
                if (centers.Count >= maxPointCount)
                    break;

                var totalMoves = Enumerable.Repeat(Vector3.Zero, centers.Count).ToList();
                var collisionCounts = Enumerable.Repeat(0.0, centers.Count).ToList();

                // ========== 创建镜像点并插入 KDTree ==========
                var kdPoints = new List<double[]>();
                var kdValues = new List<int>();
                var mirroredPoints = new List<Vector3>();
                var originalIndices = new List<int>();

                for (int dx = -1; dx <= 1; dx++)
                {
                    for (int dy = -1; dy <= 1; dy++)
                    {
                        for (int i = 0; i < centers.Count; i++)
                        {
                            var pt = new Vector3(
                                (float)(centers[i].X + dx * _tileWidth),
                                (float)(centers[i].Y + dy * _tileHeight),
                                centers[i].Z);

                            kdPoints.Add(new double[] { pt.X, pt.Y });
                            kdValues.Add(mirroredPoints.Count);
                            mirroredPoints.Add(pt);
                            originalIndices.Add(i);
                        }
                    }
                }

                var tree = new KDTree<double, int>(
                    2,
                    kdPoints.ToArray(),
                    kdValues.ToArray(),
                    (a, b) =>
                    {
                        double dx = a[0] - b[0];
                        double dy = a[1] - b[1];
                        return Math.Sqrt(dx * dx + dy * dy);
                    }
                );

                // ========== 斥力计算（保持与原逻辑一致） ==========
                for (int i = 0; i < centers.Count; i++)
                {
                    double maxSearchRadius = baseDist * _maxFactor;
                    var neighbors = tree.RadialSearch(
                        new double[] { centers[i].X, centers[i].Y },
                        maxSearchRadius);

                    foreach (var (point, idx) in neighbors)
                    {
                        int j = originalIndices[idx];
                        if (j == i) continue; // 不与自身的镜像作用

                        var mirrorJ = mirroredPoints[idx];
                        var delta = centers[i] - mirrorJ;
                        double d = delta.Length();

                        if (d < 0.001) continue;

                        double factorI = EvaluateDensityFactor(centers[i], repellers, repellerFactors);
                        double factorJ = EvaluateDensityFactor(centers[j], repellers, repellerFactors);
                        double localDist = baseDist * 0.5 * (factorI + factorJ);
                        if (d > localDist) continue;

                        double pushStrength = 0.5 * (localDist - d);
                        double maxStep = baseDist * 0.5;
                        pushStrength = Math.Min(pushStrength, maxStep);

                        var move = Vector3.Normalize(delta) * (float)pushStrength;
                        totalMoves[i] += move;
                        totalMoves[j] -= move;
                        collisionCounts[i] += 1.0;
                        collisionCounts[j] += 1.0;
                    }
                }

                for (int i = 0; i < centers.Count; i++)
                {
                    if (collisionCounts[i] > 0.0)
                        centers[i] += totalMoves[i] / (float)collisionCounts[i];
                }

                // ========== 插值生长（覆盖首尾段，按降序批量插入） ==========
                if (centers.Count < maxPointCount)
                {
                    int n = centers.Count;
                    var toInsert = new List<(int insertAt, Vector3 pt)>();

                    foreach (var (a, b) in TopologyHelpers.EnumerateSegments(n, isClosed))
                    {
                        double factorA = EvaluateDensityFactor(centers[a], repellers, repellerFactors);
                        double factorB = EvaluateDensityFactor(centers[b], repellers, repellerFactors);
                        double insertThreshold = baseDist * 0.5 * (factorA + factorB) - 1.0;

                        if (Vector3.Distance(centers[a], centers[b]) > insertThreshold)
                        {
                            var mid = 0.5f * (centers[a] + centers[b]);
                            // 约定将新点插在段 (a→b) 的 b 索引处
                            toInsert.Add((b, mid));
                        }
                    }

                    if (toInsert.Count > 0)
                    {
                        // 关键：按索引降序插入，避免因插入导致后续索引位移
                        toInsert.Sort((x, y) => y.insertAt.CompareTo(x.insertAt));

                        foreach (var (insertAt, pt) in toInsert)
                        {
                            if (centers.Count >= maxPointCount) break;
                            centers.Insert(insertAt, pt);
                        }
                    }
                }
            }

            return centers;
        }

        private double EvaluateDensityFactor(Vector3 pt, List<Vector3> repellers, List<double> factors)
        {
            if (repellers.Count == 0) return 1.0;
            double maxLocalFactor = 1.0;

            for (int i = 0; i < repellers.Count; i++)
            {
                double d = Vector3.Distance(pt, repellers[i]);
                if (d > _maxEffectDist) continue;

                double t = d / _maxEffectDist;
                double repellerStrength = factors[Math.Min(i, factors.Count - 1)];
                double localFactor = 1.0 + (_maxFactor - 1.0) * repellerStrength * (1.0 - t);

                if (localFactor > maxLocalFactor)
                    maxLocalFactor = localFactor;
            }

            return maxLocalFactor;
        }
    }
}
