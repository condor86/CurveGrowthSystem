using Supercluster.KDTree;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

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
                
                if (iter == 0)
                {
                    // 根目录路径（即可执行文件所在目录）
                    string rootDir = AppDomain.CurrentDomain.BaseDirectory;

                    // 上级目录，用于保存输出结果
                    string parentDir = Path.GetFullPath(Path.Combine(rootDir, "..", "..", ".."));

                    string mirroredPath = Path.Combine(parentDir, "resultsMirroredFlat.csv");
                    string indexPath = Path.Combine(parentDir, "resultsOriginalIndices.csv");

                    // 输出 mirroredFlat 为 CSV，每行为 x,y,z
                    using (var writer = new StreamWriter(mirroredPath))
                    {
                        for (int i = 0; i < mirroredPoints.Count; i++)
                        {
                            var pt = mirroredPoints[i];
                            writer.WriteLine($"{pt.X},{pt.Y},{pt.Z}");
                        }
                    }

                    // 输出 originalIndices 为 CSV，每行一个整数
                    using (var writer = new StreamWriter(indexPath))
                    {
                        for (int i = 0; i < originalIndices.Count; i++)
                        {
                            writer.WriteLine(originalIndices[i]);
                        }
                    }

                    Console.WriteLine($"初始镜像点与索引已输出：\n→ {mirroredPath}\n→ {indexPath}");
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
                
                // === 满足终止条件时输出 totalMoves 与 collisionCounts ===
                if (iter == 1)
                {
                    // 根目录路径（即可执行文件所在目录）
                    string rootDir = AppDomain.CurrentDomain.BaseDirectory;

                    // 上级目录，用于保存输出结果
                    string parentDir = Path.GetFullPath(Path.Combine(rootDir, "..", "..", ".."));
    
                    string movePath = Path.Combine(parentDir, "resultsTotalMoves.csv");
                    string countPath = Path.Combine(parentDir, "resultsCollisionCounts.csv");

                    // 输出 totalMoves 为 CSV，每行为 x,y,z
                    using (var writer = new StreamWriter(movePath))
                    {
                        for (int i = 0; i < totalMoves.Count; i++)
                        {
                            Vector3 move = totalMoves[i];
                            writer.WriteLine($"{move.X},{move.Y},{move.Z}");
                        }
                    }

                    // 输出 collisionCounts 为 CSV，每行一个值
                    using (var writer = new StreamWriter(countPath))
                    {
                        for (int i = 0; i < collisionCounts.Count; i++)
                        {
                            double val = collisionCounts[i];
                            writer.WriteLine($"{val}");
                        }
                    }

                    Console.WriteLine($"测试结果已保存至：{movePath}");
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
