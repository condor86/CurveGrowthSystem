using Supercluster.KDTree;
using System;
using System.Collections.Generic;
using NumSharp;

namespace CrvGrowth
{
    public class GrowthSystem
    {
        private readonly double _tileWidth = 1000.0;
        private readonly double _tileHeight = 1000.0;
        private readonly double _maxFactor = 1.5;
        private readonly double _maxEffectDist = 300.0;

        public NDArray Run(
            NDArray starting,           // [N, 3]
            NDArray repellers,          // [M, 3]
            NDArray repellerFactors,    // [M]
            int maxPointCount = 200,
            int maxIterCount = 200,
            double baseDist = 75.0)
        {
            var centers = starting.copy();  // [N, 3]

            for (int iter = 0; iter < maxIterCount; iter++)
            {
                int N = centers.shape[0];
                if (N >= maxPointCount) break;

                var totalMoves = np.zeros_like(centers);  // [N, 3]
                var collisionCounts = np.zeros(N);        // [N]

                // === 构造镜像点 ===
                var offsets = np.array(new double[,]
                {
                    {-_tileWidth, -_tileHeight},
                    { 0,         -_tileHeight},
                    {_tileWidth, -_tileHeight},
                    {-_tileWidth, 0},
                    { 0,          0},
                    {_tileWidth,  0},
                    {-_tileWidth, _tileHeight},
                    { 0,          _tileHeight},
                    {_tileWidth,  _tileHeight}
                });  // [9, 2]

                var offset9N3 = np.zeros((9, N, 3));
                offset9N3[":", ":", 0] = offsets[":", 0].reshape(9, 1);
                offset9N3[":", ":", 1] = offsets[":", 1].reshape(9, 1);

                var centers9N3 = np.zeros(new Shape(9, N, 3));
                for (int i = 0; i < 9; i++)
                    centers9N3[i, ":", ":"] = centers[":", ":"];

                var mirroredFlat = (centers9N3 + offset9N3).reshape(9 * N, 3); // [9N, 3]

                // === 构造 KDTree 数据 ===
                var kdPoints = new List<double[]>();
                var kdValues = new List<int>();
                var originalIndices = new List<int>();

                for (int i = 0; i < mirroredFlat.shape[0]; i++)
                {
                    double x = (double)mirroredFlat[i, 0];
                    double y = (double)mirroredFlat[i, 1];
                    kdPoints.Add(new double[] { x, y });
                    kdValues.Add(i);
                    originalIndices.Add(i % N);
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
                    });

                // === 斥力运算 ===
                for (int i = 0; i < N; i++)
                {
                    var centerPt = centers[i, ":"];
                    var neighbors = tree.RadialSearch(
                        new double[] {
                            (double)centerPt[0],
                            (double)centerPt[1]
                        },
                        baseDist * _maxFactor);

                    foreach (var (pt, idx) in neighbors)
                    {
                        int j = originalIndices[idx];
                        if (j == i) continue;

                        var mirrorJ = mirroredFlat[idx, ":"];
                        var delta = centerPt - mirrorJ;

                        double d = Math.Sqrt(
                            Math.Pow((double)delta[0], 2) +
                            Math.Pow((double)delta[1], 2) +
                            Math.Pow((double)delta[2], 2)
                        );
                        if (d < 0.001) continue;

                        double factorI = EvaluateDensityFactor(centerPt, repellers, repellerFactors);
                        double factorJ = EvaluateDensityFactor(centers[j, ":"], repellers, repellerFactors);
                        double localDist = baseDist * 0.5 * (factorI + factorJ);
                        if (d > localDist) continue;

                        double pushStrength = 0.5 * (localDist - d);
                        double maxStep = baseDist * 0.5;
                        pushStrength = Math.Min(pushStrength, maxStep);

                        var move = delta / d * pushStrength;
                        totalMoves[i, ":"] += move;
                        totalMoves[j, ":"] -= move;
                        collisionCounts[i] += 1.0;
                        collisionCounts[j] += 1.0;
                    }
                }

                // === 应用运动 ===
                for (int i = 0; i < N; i++)
                {
                    if (collisionCounts[i] > 0.0)
                        centers[i, ":"] += totalMoves[i, ":"] / collisionCounts[i];
                }

                // === 插值添加点 ===
                N = centers.shape[0];
                if (N >= maxPointCount) break;

                var splitIndices = new List<int>();
                for (int i = 0; i < N - 1; i++)
                {
                    double factorA = EvaluateDensityFactor(centers[i, ":"], repellers, repellerFactors);
                    double factorB = EvaluateDensityFactor(centers[i + 1, ":"], repellers, repellerFactors);
                    double insertThreshold = baseDist * 0.5 * (factorA + factorB) - 1.0;

                    var diff = centers[i + 1, ":"] - centers[i, ":"];
                    double d = Math.Sqrt(
                        Math.Pow((double)diff[0], 2) +
                        Math.Pow((double)diff[1], 2) +
                        Math.Pow((double)diff[2], 2)
                    );

                    if (d > insertThreshold)
                        splitIndices.Add(i + 1 + splitIndices.Count);
                }

                if (splitIndices.Count > 0)
                {
                    var newCenters = new List<NDArray>();
                    for (int i = 0; i < N; i++)
                        newCenters.Add(centers[i, ":"].reshape(1, 3));

                    foreach (int idx in splitIndices)
                    {
                        if (newCenters.Count >= maxPointCount) break;
                        var mid = 0.5 * (newCenters[idx - 1] + newCenters[idx]);
                        newCenters.Insert(idx, mid);
                    }

                    centers = np.concatenate(newCenters.ToArray(), axis: 0);
                }
            }

            return centers;  // [N_final, 3]
        }

        private double EvaluateDensityFactor(NDArray pt, NDArray repellers, NDArray factors)
        {
            if (repellers.shape[0] == 0) return 1.0;

            double maxLocalFactor = 1.0;
            int M = repellers.shape[0];

            for (int i = 0; i < M; i++)
            {
                double dx = (double)repellers[i, 0] - (double)pt[0];
                double dy = (double)repellers[i, 1] - (double)pt[1];
                double dz = (double)repellers[i, 2] - (double)pt[2];
                double dist = Math.Sqrt(dx * dx + dy * dy + dz * dz);

                if (dist > _maxEffectDist) continue;

                double t = dist / _maxEffectDist;
                double strength = (double)factors[i];
                double localFactor = 1.0 + (_maxFactor - 1.0) * strength * (1.0 - t);

                if (localFactor > maxLocalFactor)
                    maxLocalFactor = localFactor;
            }

            return maxLocalFactor;
        }
    }
}
