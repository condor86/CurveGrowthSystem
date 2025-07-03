using Supercluster.KDTree;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using NumSharp;

namespace CrvGrowth
{
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
            double baseDist = 75.0)
        {
            var centers = new List<Vector3>(starting);

            for (int iter = 0; iter < maxIterCount; iter++)
            {
                if (centers.Count >= maxPointCount)
                    break;

                var totalMoves = Enumerable.Repeat(Vector3.Zero, centers.Count).ToList();
                var collisionCounts = Enumerable.Repeat(0.0, centers.Count).ToList();

                // ========== 使用 NumSharp 向量化创建镜像点 ==========
                int N = centers.Count;
                double[,] centerArray = new double[N, 3];
                for (int i = 0; i < N; i++)
                {
                    centerArray[i, 0] = centers[i].X;
                    centerArray[i, 1] = centers[i].Y;
                    centerArray[i, 2] = centers[i].Z;
                }
                NDArray centerND = np.array(centerArray);  // [N, 3]
//                Console.WriteLine("centerND:");
//                Console.WriteLine(centerND.ToString(true));
                
                NDArray offsets = np.array(new double[,]
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

                NDArray offset9N3 = np.zeros((9, N, 3));
                offset9N3[$":", $":", 0] = offsets[$":", 0].reshape(9, 1);
                offset9N3[$":", $":", 1] = offsets[$":", 1].reshape(9, 1);
                
//                Console.WriteLine("offsets:\n" + offsets.ToString());
//                Console.WriteLine("offset9N3:\n" + offset9N3.ToString());
                // 创建 shape 为 [9, N, 3] 的空数组
                NDArray centers9N3 = np.zeros((9, N, 3));

                for (int i = 0; i < 9; i++)  // 9 个镜像方向
                {
                    for (int j = 0; j < N; j++)  // N 个点
                    {
                        centers9N3[i, j, 0] = centerND[j, 0];  // X
                        centers9N3[i, j, 1] = centerND[j, 1];  // Y
                        centers9N3[i, j, 2] = centerND[j, 2];  // Z
                    }
                }
//                NDArray centers9N3 = centerND;
//                NDArray centers9N3 = np.tile(centerND, (9, 1)).reshape(9, N, 3);
                NDArray mirroredFlat = (centers9N3 + offset9N3).reshape(N * 9, 3);

                var kdPoints = new List<double[]>();
                var mirroredPoints = new List<Vector3>();
                var kdValues = new List<int>();
                var originalIndices = new List<int>();

                for (int i = 0; i < mirroredFlat.shape[0]; i++)
                {
                    double x = (double)mirroredFlat[i, 0];
                    double y = (double)mirroredFlat[i, 1];
                    double z = (double)mirroredFlat[i, 2];

                    kdPoints.Add(new double[] { x, y });
                    mirroredPoints.Add(new Vector3((float)x, (float)y, (float)z));
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
                    }
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

                // ========== 插值生长 ==========
                if (centers.Count < maxPointCount)
                {
                    var splitIndices = new List<int>();

                    for (int i = 0; i < centers.Count - 1; i++)
                    {
                        double factorA = EvaluateDensityFactor(centers[i], repellers, repellerFactors);
                        double factorB = EvaluateDensityFactor(centers[i + 1], repellers, repellerFactors);
                        double insertThreshold = baseDist * 0.5 * (factorA + factorB) - 1.0;

                        if (Vector3.Distance(centers[i], centers[i + 1]) > insertThreshold)
                            splitIndices.Add(i + 1 + splitIndices.Count);
                    }

                    foreach (int splitIndex in splitIndices)
                    {
                        var a = centers[splitIndex - 1];
                        var b = centers[splitIndex];
                        var newCenter = 0.5f * (a + b);
                        centers.Insert(splitIndex, newCenter);

                        if (centers.Count >= maxPointCount) break;
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
