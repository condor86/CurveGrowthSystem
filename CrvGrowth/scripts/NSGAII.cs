// File: NSGAII.cs
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace NSGAII
{
    /// <summary>
    /// 配置参数（与 Wallacei 对齐：PopulationSize, Generations, CrossoverRate, MutationRate, η 等）
    /// 必填：
    ///  - GeneLength
    ///  - LowerBounds / UpperBounds（长度 == GeneLength）
    ///  - Evaluate（Func<double[], double[]>，返回 {f1, f2, ...}，默认“越小越好”）
    /// </summary>
    public class NSGAConfig
    {
        public int    PopulationSize        = 50;
        public int    Generations           = 100;
        public double CrossoverRate         = 0.9;         // 交叉概率
        public double MutationRate          = 0.05;        // 变异概率（常用 1/n）
        public int    GeneLength            = 0;           // 必填
        public double[] LowerBounds         = null!;       // 必填
        public double[] UpperBounds         = null!;       // 必填

        public int    RandomSeed            = 1;
        public int?   DegreeOfParallelism   = null;        // Evaluate 并行评估
        public string? LogDir               = null;        // 每代导出 front0 / bestGenes.csv

        // SBX 与 Polynomial 变异的分布指数（distribution index）
        public double SbxEta                = 20.0;        // η_c
        public double PolyMutationEta       = 20.0;        // η_m

        /// <summary>
        /// 评价函数：输入 genes，返回目标值数组（默认“越小越好”）。
        /// 如有“越大越好”的目标，请在此处取负以统一为最小化。
        /// </summary>
        public Func<double[], double[]> Evaluate = null!;
    }

    /// <summary>
    /// 个体
    /// </summary>
    public class Individual
    {
        public double[] Genes;
        public double[] Objectives;
        public int Rank;
        public double CrowdingDistance;

        public Individual(int geneLength, int nObjectives = 2)
        {
            Genes = new double[geneLength];
            Objectives = new double[nObjectives];
            Rank = int.MaxValue;
            CrowdingDistance = 0.0;
        }

        public Individual Clone()
        {
            var c = new Individual(Genes.Length, Objectives.Length);
            Genes.CopyTo(c.Genes, 0);
            Objectives.CopyTo(c.Objectives, 0);
            c.Rank = Rank;
            c.CrowdingDistance = CrowdingDistance;
            return c;
        }
    }

    /// <summary>
    /// NSGA-II 主类：初始化→评估→SBX交叉→Polynomial变异→合并→非支配排序→拥挤度→选择→循环
    /// </summary>
    public class NSGAII
    {
        private readonly NSGAConfig _cfg;
        private readonly Random _rand;

        public NSGAII(NSGAConfig cfg)
        {
            _cfg = cfg ?? throw new ArgumentNullException(nameof(cfg));

            if (_cfg.GeneLength <= 0)
                throw new ArgumentException("GeneLength must be > 0");
            if (_cfg.LowerBounds == null || _cfg.UpperBounds == null)
                throw new ArgumentException("LowerBounds/UpperBounds must be provided.");
            if (_cfg.LowerBounds.Length != _cfg.GeneLength || _cfg.UpperBounds.Length != _cfg.GeneLength)
                throw new ArgumentException("Bounds length must equal GeneLength.");
            if (_cfg.Evaluate == null)
                throw new ArgumentException("Evaluate callback must be provided.");

            _rand = new Random(_cfg.RandomSeed);

            if (!string.IsNullOrEmpty(_cfg.LogDir))
                Directory.CreateDirectory(_cfg.LogDir);
        }

        /// <summary>
        /// 运行 NSGA-II，返回最终一代种群（可取 Rank==0 作为 Pareto 前沿）
        /// </summary>
        public List<Individual> Run()
        {
            // 初始化
            var pop = InitializePopulation();

            // 主循环
            for (int gen = 0; gen < _cfg.Generations; gen++)
            {
                EvaluatePopulation(pop);

                var offspring = Reproduce(pop);
                EvaluatePopulation(offspring);

                var combined = pop.Concat(offspring).ToList();
                var fronts = FastNonDominatedSort(combined);

                pop = SelectNewPopulation(fronts, _cfg.PopulationSize);

                Console.WriteLine($"Gen {gen + 1}/{_cfg.Generations} | Front0 size = {fronts[0].Count}");

                // 日志
                if (!string.IsNullOrEmpty(_cfg.LogDir))
                {
                    ExportFrontToCsv(fronts[0], Path.Combine(_cfg.LogDir!, $"gen_{gen + 1}_front0.csv"));
                    ExportBestGenes(fronts[0], Path.Combine(_cfg.LogDir!, $"gen_{gen + 1}_bestGenes.csv"));
                }
            }

            return pop;
        }

        // ============================ 初始化 / 评估 ============================

        private List<Individual> InitializePopulation()
        {
            var list = new List<Individual>(_cfg.PopulationSize);
            for (int i = 0; i < _cfg.PopulationSize; i++)
            {
                var ind = new Individual(_cfg.GeneLength, nObjectives: 2); // 默认双目标；如需更多，可在 Evaluate 返回更多维度
                for (int g = 0; g < _cfg.GeneLength; g++)
                {
                    double lo = _cfg.LowerBounds[g];
                    double hi = _cfg.UpperBounds[g];
                    ind.Genes[g] = lo + _rand.NextDouble() * (hi - lo);
                }
                list.Add(ind);
            }
            return list;
        }

        private void EvaluatePopulation(List<Individual> pop)
        {
            if (_cfg.DegreeOfParallelism.HasValue && _cfg.DegreeOfParallelism.Value > 1)
            {
                var opts = new ParallelOptions { MaxDegreeOfParallelism = _cfg.DegreeOfParallelism.Value };
                Parallel.ForEach(pop, opts, ind =>
                {
                    var f = _cfg.Evaluate(ind.Genes);
                    EnsureObjectiveSize(ind, f.Length);
                    Buffer.BlockCopy(f, 0, ind.Objectives, 0, sizeof(double) * f.Length);
                });
            }
            else
            {
                foreach (var ind in pop)
                {
                    var f = _cfg.Evaluate(ind.Genes);
                    EnsureObjectiveSize(ind, f.Length);
                    Buffer.BlockCopy(f, 0, ind.Objectives, 0, sizeof(double) * f.Length);
                }
            }
        }

        private static void EnsureObjectiveSize(Individual ind, int n)
        {
            if (ind.Objectives.Length != n)
                ind.Objectives = new double[n];
        }

        // ============================ 遗传操作 ============================

        private List<Individual> Reproduce(List<Individual> pop)
        {
            var children = new List<Individual>(_cfg.PopulationSize);
            while (children.Count < _cfg.PopulationSize)
            {
                var p1 = TournamentSelect(pop);
                var p2 = TournamentSelect(pop);

                Individual c1, c2;
                CrossoverSBX(p1, p2, out c1, out c2);
                MutatePolynomial(c1);
                MutatePolynomial(c2);
                ApplyBounds(c1);
                ApplyBounds(c2);

                children.Add(c1);
                if (children.Count < _cfg.PopulationSize) children.Add(c2);
            }
            return children;
        }

        private Individual TournamentSelect(List<Individual> pop)
        {
            var a = pop[_rand.Next(pop.Count)];
            var b = pop[_rand.Next(pop.Count)];
            if (a.Rank < b.Rank) return a;
            if (a.Rank > b.Rank) return b;
            return a.CrowdingDistance >= b.CrowdingDistance ? a : b;
        }

        /// <summary>
        /// SBX（Simulated Binary Crossover），按维度独立；参考 Deb (2002)
        /// </summary>
        private void CrossoverSBX(Individual p1, Individual p2, out Individual c1, out Individual c2)
        {
            c1 = p1.Clone();
            c2 = p2.Clone();
            if (_rand.NextDouble() >= _cfg.CrossoverRate) return;

            for (int i = 0; i < _cfg.GeneLength; i++)
            {
                double x1 = p1.Genes[i];
                double x2 = p2.Genes[i];

                if (Math.Abs(x1 - x2) < 1e-14)
                {
                    c1.Genes[i] = x1;
                    c2.Genes[i] = x2;
                    continue;
                }

                double y1 = Math.Min(x1, x2);
                double y2 = Math.Max(x1, x2);
                double lo = _cfg.LowerBounds[i];
                double hi = _cfg.UpperBounds[i];

                double rand = _rand.NextDouble();
                double beta = 1.0 + (2.0 * (y1 - lo) / (y2 - y1));
                double alpha = 2.0 - Math.Pow(beta, -(_cfg.SbxEta + 1.0));
                double betaq;

                if (rand <= 1.0 / alpha)
                    betaq = Math.Pow(rand * alpha, 1.0 / (_cfg.SbxEta + 1.0));
                else
                    betaq = Math.Pow(1.0 / (2.0 - rand * alpha), 1.0 / (_cfg.SbxEta + 1.0));

                double child1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
                double child2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));

                // 随机交换左右
                if (_rand.NextDouble() < 0.5)
                {
                    c1.Genes[i] = Clamp(child1, lo, hi);
                    c2.Genes[i] = Clamp(child2, lo, hi);
                }
                else
                {
                    c1.Genes[i] = Clamp(child2, lo, hi);
                    c2.Genes[i] = Clamp(child1, lo, hi);
                }
            }
        }

        /// <summary>
        /// 多项式变异（Polynomial Mutation）；参考 Deb (2002)
        /// </summary>
        private void MutatePolynomial(Individual ind)
        {
            for (int i = 0; i < _cfg.GeneLength; i++)
            {
                if (_rand.NextDouble() >= _cfg.MutationRate) continue;

                double lo = _cfg.LowerBounds[i];
                double hi = _cfg.UpperBounds[i];
                double x = ind.Genes[i];

                if (hi - lo < 1e-14) continue;

                double u = _rand.NextDouble();
                double delta1 = (x - lo) / (hi - lo);
                double delta2 = (hi - x) / (hi - lo);
                double mutPow = 1.0 / (_cfg.PolyMutationEta + 1.0);
                double deltaq;

                if (u < 0.5)
                {
                    double xy = 1.0 - delta1;
                    double val = 2.0 * u + (1.0 - 2.0 * u) * Math.Pow(xy, _cfg.PolyMutationEta + 1.0);
                    deltaq = Math.Pow(val, mutPow) - 1.0;
                }
                else
                {
                    double xy = 1.0 - delta2;
                    double val = 2.0 * (1.0 - u) + 2.0 * (u - 0.5) * Math.Pow(xy, _cfg.PolyMutationEta + 1.0);
                    deltaq = 1.0 - Math.Pow(val, mutPow);
                }

                x = x + deltaq * (hi - lo);
                ind.Genes[i] = Clamp(x, lo, hi);
            }
        }

        private void ApplyBounds(Individual ind)
        {
            for (int i = 0; i < _cfg.GeneLength; i++)
            {
                if (ind.Genes[i] < _cfg.LowerBounds[i]) ind.Genes[i] = _cfg.LowerBounds[i];
                if (ind.Genes[i] > _cfg.UpperBounds[i]) ind.Genes[i] = _cfg.UpperBounds[i];
            }
        }

        private static double Clamp(double v, double lo, double hi)
            => v < lo ? lo : (v > hi ? hi : v);

        // ============================ 非支配排序 / 拥挤度 / 选择 ============================

        /// <summary>
        /// 快速非支配排序
        /// </summary>
        private List<List<Individual>> FastNonDominatedSort(List<Individual> pop)
        {
            var fronts = new List<List<Individual>>();
            var S = new Dictionary<Individual, List<Individual>>(pop.Count);
            var n = new Dictionary<Individual, int>(pop.Count);

            var F1 = new List<Individual>();

            foreach (var p in pop)
            {
                S[p] = new List<Individual>();
                n[p] = 0;

                foreach (var q in pop)
                {
                    if (ReferenceEquals(p, q)) continue;

                    if (Dominates(p, q))
                        S[p].Add(q);
                    else if (Dominates(q, p))
                        n[p]++;
                }

                if (n[p] == 0)
                {
                    p.Rank = 0;
                    F1.Add(p);
                }
            }

            fronts.Add(F1);
            int i = 0;
            while (i < fronts.Count && fronts[i].Count > 0)
            {
                var next = new List<Individual>();
                foreach (var p in fronts[i])
                {
                    foreach (var q in S[p])
                    {
                        n[q]--;
                        if (n[q] == 0)
                        {
                            q.Rank = i + 1;
                            next.Add(q);
                        }
                    }
                }
                i++;
                fronts.Add(next);
            }

            return fronts;
        }

        /// <summary>
        /// p 是否支配 q（默认所有目标“越小越好”）
        /// </summary>
        private static bool Dominates(Individual p, Individual q)
        {
            bool betterInAny = false;
            int m = p.Objectives.Length;

            for (int i = 0; i < m; i++)
            {
                if (p.Objectives[i] > q.Objectives[i]) return false; // 任一目标上更差则不支配
                if (p.Objectives[i] < q.Objectives[i]) betterInAny = true;
            }
            return betterInAny;
        }

        /// <summary>
        /// 拥挤距离计算
        /// </summary>
        private static void CalculateCrowdingDistance(List<Individual> front)
        {
            if (front.Count == 0) return;

            foreach (var ind in front)
                ind.CrowdingDistance = 0.0;

            int m = front[0].Objectives.Length;

            for (int obj = 0; obj < m; obj++)
            {
                var sorted = front.OrderBy(ind => ind.Objectives[obj]).ToList();

                // 边界点设为无穷大，优先保留
                sorted[0].CrowdingDistance = double.PositiveInfinity;
                sorted[^1].CrowdingDistance = double.PositiveInfinity;

                double min = sorted[0].Objectives[obj];
                double max = sorted[^1].Objectives[obj];
                double range = max - min;
                if (range <= 1e-14) continue;

                for (int i = 1; i < sorted.Count - 1; i++)
                {
                    double next = sorted[i + 1].Objectives[obj];
                    double prev = sorted[i - 1].Objectives[obj];
                    sorted[i].CrowdingDistance += (next - prev) / range;
                }
            }
        }

        /// <summary>
        /// 环境选择：按前沿依次装入，最后一层按拥挤距离截断
        /// </summary>
        private static List<Individual> SelectNewPopulation(List<List<Individual>> fronts, int targetSize)
        {
            var newPop = new List<Individual>(targetSize);
            int i = 0;

            while (i < fronts.Count && newPop.Count + fronts[i].Count <= targetSize)
            {
                CalculateCrowdingDistance(fronts[i]);
                newPop.AddRange(fronts[i]);
                i++;
            }

            if (newPop.Count < targetSize && i < fronts.Count && fronts[i].Count > 0)
            {
                CalculateCrowdingDistance(fronts[i]);
                var rest = fronts[i].OrderByDescending(x => x.CrowdingDistance)
                                    .Take(targetSize - newPop.Count);
                newPop.AddRange(rest);
            }

            return newPop;
        }

        // ============================ 日志导出 ============================

        private void ExportFrontToCsv(List<Individual> front, string path)
        {
            using var sw = new StreamWriter(path, false);
            // header: f0,f1,...,g0,g1,...
            var fCols = Enumerable.Range(0, front[0].Objectives.Length).Select(i => $"f{i}");
            var gCols = Enumerable.Range(0, _cfg.GeneLength).Select(i => $"g{i}");
            sw.WriteLine(string.Join(",", fCols.Concat(gCols)));

            foreach (var ind in front)
            {
                var f = string.Join(",", ind.Objectives.Select(v => v.ToString("G17")));
                var g = string.Join(",", ind.Genes.Select(v => v.ToString("G17")));
                sw.WriteLine($"{f},{g}");
            }
        }

        /// <summary>
        /// 以 f 向量的 L1 和为排序指标，导出一个代表解的 genes（仅日志用途）
        /// </summary>
        private void ExportBestGenes(List<Individual> front, string path)
        {
            var best = front.OrderBy(ind => ind.Objectives.Sum()).First();
            using var sw = new StreamWriter(path, false);
            sw.WriteLine(string.Join(",", best.Genes.Select(v => v.ToString("G17"))));
        }
    }
}
