using Supercluster.KDTree;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CrvGrowth
{
    public class GrowthSystem
    {
        private readonly double _tileWidth = 1000.0;
        private readonly double _tileHeight = 1000.0;
        private readonly double _maxFactor = 1.5;
        private readonly double _maxEffectDist = 300.0;

        public List<Point3D> Run(
            List<Point3D> starting,
            List<Point3D> repellers,
            List<double> repellerFactors,
            int maxPointCount = 200,
            int maxIterCount = 200,
            double baseDist = 75.0)
        {
            var centers = new List<Point3D>(starting);

            for (int iter = 0; iter < maxIterCount; iter++)
            {
                if (centers.Count >= maxPointCount)
                    break;

                var totalMoves = Enumerable.Repeat(Vector3D.Zero, centers.Count).ToList();
                var collisionCounts = Enumerable.Repeat(0.0, centers.Count).ToList();

                // ========== 创建镜像点并插入 KDTree ==========
                var kdPoints = new List<double[]>();
                var kdValues = new List<int>();
                var mirroredPoints = new List<Point3D>();
                var originalIndices = new List<int>();

                for (int dx = -1; dx <= 1; dx++)
                {
                    for (int dy = -1; dy <= 1; dy++)
                    {
                        for (int i = 0; i < centers.Count; i++)
                        {
                            var pt = new Point3D(
                                centers[i].X + dx * _tileWidth,
                                centers[i].Y + dy * _tileHeight,
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
                    new Func<double[], double[], double>(delegate (double[] a, double[] b)
                    {
                        double dx = a[0] - b[0];
                        double dy = a[1] - b[1];
                        return Math.Sqrt(dx * dx + dy * dy);
                    })
                );

                // ========== 斥力计算 ==========
                for (int i = 0; i < centers.Count; i++)
                {
                    double maxSearchRadius = baseDist * _maxFactor;
                    var neighbors = tree.RadialSearch(
                        new double[] { centers[i].X, centers[i].Y },
                        maxSearchRadius);

                    foreach (var (point, idx) in neighbors)
                    {
                        int j = originalIndices[idx];
                        if (j == i) continue;

                        var mirrorJ = mirroredPoints[idx];
                        var delta = centers[i] - mirrorJ;
                        double d = delta.Length;
                        if (d < 0.001) continue;

                        double factorI = EvaluateDensityFactor(centers[i], repellers, repellerFactors);
                        double factorJ = EvaluateDensityFactor(centers[j], repellers, repellerFactors);
                        double localDist = baseDist * 0.5 * (factorI + factorJ);
                        if (d > localDist) continue;

                        double pushStrength = 0.5 * (localDist - d);
                        double maxStep = baseDist * 0.5;
                        pushStrength = Math.Min(pushStrength, maxStep);

                        var move = delta.Normalized * pushStrength;
                        totalMoves[i] += move;
                        totalMoves[j] -= move;
                        collisionCounts[i] += 1.0;
                        collisionCounts[j] += 1.0;
                    }
                }

                for (int i = 0; i < centers.Count; i++)
                {
                    if (collisionCounts[i] > 0.0)
                        centers[i] += totalMoves[i] / collisionCounts[i];
                }

                // ========== 插值生长 ==========
                if (centers.Count < maxPointCount)
                {
                    var splitIndices = new List<int>();

                    for (int i = 0; i < centers.Count - 1; i++)
                    {
                        double factorA = EvaluateDensityFactor(centers[i], repellers, repellerFactors);
                        double factorB = EvaluateDensityFactor(centers[i + 1], repellers, repellerFactors);
                        double insertThreshold = baseDist * 0.5 * (factorA + factorB) - 1.0;

                        if ((centers[i] - centers[i + 1]).Length > insertThreshold)
                            splitIndices.Add(i + 1 + splitIndices.Count);
                    }

                    foreach (int splitIndex in splitIndices)
                    {
                        var a = centers[splitIndex - 1];
                        var b = centers[splitIndex];
                        var newCenter = new Point3D(
                            0.5 * (a.X + b.X),
                            0.5 * (a.Y + b.Y),
                            0.5 * (a.Z + b.Z));
                        centers.Insert(splitIndex, newCenter);

                        if (centers.Count >= maxPointCount) break;
                    }
                }
            }

            return centers;
        }

        private double EvaluateDensityFactor(Point3D pt, List<Point3D> repellers, List<double> factors)
        {
            if (repellers.Count == 0) return 1.0;
            double maxLocalFactor = 1.0;

            for (int i = 0; i < repellers.Count; i++)
            {
                double d = (pt - repellers[i]).Length;
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
