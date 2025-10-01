// File: NSGAII.cs
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Diagnostics;      // ★ 新增：用于计时
using CrvGrowth;               // 使用 NSGAEvalContext（从 NSGAWiring.cs 提供）

namespace NSGAII
{
    public class NSGAConfig
    {
        public int    PopulationSize        = 50;
        public int    Generations           = 100;
        public double CrossoverRate         = 0.9;
        public double MutationRate          = 0.05;        // 常用 1.0 / GeneLength
        public int    GeneLength            = 0;           // 必填
        public double[] LowerBounds         = null!;
        public double[] UpperBounds         = null!;

        public int    RandomSeed            = 1;
        public int?   DegreeOfParallelism   = null;        // 并行评估的最大并行度
        public string? LogDir               = null;        // 每代导出 front0 / bestGenes.csv

        public double SbxEta                = 20.0;        // SBX η_c
        public double PolyMutationEta       = 20.0;        // 变异 η_m

        /// <summary>评价函数：返回目标数组（默认“越小越好”）。若某目标是“越大越好”，请在外部取负。</summary>
        public Func<double[], double[]> Evaluate = null!;
    }

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

        public List<Individual> Run()
        {
            var wall = Stopwatch.StartNew();   // ★ 总耗时计时器

            // 初始化
            var pop = InitializePopulation();

            // 主循环
            for (int gen = 0; gen < _cfg.Generations; gen++)
            {
                int genHuman = gen + 1; // 1-based 代数

                // 评估父代
                EvaluatePopulation(pop, genHuman);

                // 生成子代并评估
                var offspring = Reproduce(pop);
                EvaluatePopulation(offspring, genHuman);

                // 合并 → 非支配排序 → 选择
                var combined = pop.Concat(offspring).ToList();
                var fronts = FastNonDominatedSort(combined);
                pop = SelectNewPopulation(fronts, _cfg.PopulationSize);

                // ★ 计算该代的“当前最优（夏最小 / 冬最大）”
                (double bestSummer, double bestWinter) = GetGenerationMetrics(pop);

                Console.WriteLine(
                    $"Gen {genHuman}/{_cfg.Generations} | Front0 size = {fronts[0].Count} | " +
                    $"Elapsed = {wall.Elapsed} | Best(Summer min, Winter max) = ({bestSummer:G6}, {bestWinter:G6})");

                // 日志导出
                if (!string.IsNullOrEmpty(_cfg.LogDir))
                {
                    ExportFrontToCsv(fronts[0], Path.Combine(_cfg.LogDir!, $"gen_{genHuman}_front0.csv"));
                    ExportBestGenes(fronts[0], Path.Combine(_cfg.LogDir!, $"gen_{genHuman}_bestGenes.csv"));
                }
            }

            wall.Stop();
            return pop;
        }

        // ===== 初始化 / 评估 =====

        private List<Individual> InitializePopulation()
        {
            var list = new List<Individual>(_cfg.PopulationSize);
            for (int i = 0; i < _cfg.PopulationSize; i++)
            {
                var ind = new Individual(_cfg.GeneLength, nObjectives: 2); // 默认双目标
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

        /// <summary>评估一个种群（设置“第几代｜第几个个体”的上下文），支持并行。</summary>
        private void EvaluatePopulation(List<Individual> pop, int genHuman)
        {
            int maxPar = _cfg.DegreeOfParallelism ?? -1;

            if (maxPar > 1 || maxPar == -1)
            {
                var opts = new ParallelOptions { MaxDegreeOfParallelism = maxPar };
                Parallel.For(0, pop.Count, opts, k =>
                {
                    try
                    {
                        NSGAEvalContext.Set(genHuman, k + 1); // 1-based 个体编号
                        var f = _cfg.Evaluate(pop[k].Genes);
                        EnsureObjectiveSize(pop[k], f.Length);
                        for (int i = 0; i < f.Length; i++) pop[k].Objectives[i] = f[i];
                    }
                    finally
                    {
                        NSGAEvalContext.Clear();
                    }
                });
            }
            else
            {
                for (int k = 0; k < pop.Count; k++)
                {
                    try
                    {
                        NSGAEvalContext.Set(genHuman, k + 1);
                        var f = _cfg.Evaluate(pop[k].Genes);
                        EnsureObjectiveSize(pop[k], f.Length);
                        for (int i = 0; i < f.Length; i++) pop[k].Objectives[i] = f[i];
                    }
                    finally
                    {
                        NSGAEvalContext.Clear();
                    }
                }
            }
        }

        private static void EnsureObjectiveSize(Individual ind, int n)
        {
            if (ind.Objectives.Length != n)
                ind.Objectives = new double[n];
        }

        // ===== 遗传操作 =====

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

        // ===== 非支配排序 / 拥挤度 / 选择 =====

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

        private static bool Dominates(Individual p, Individual q)
        {
            bool betterInAny = false;
            int m = p.Objectives.Length;

            for (int i = 0; i < m; i++)
            {
                if (p.Objectives[i] > q.Objectives[i]) return false; // 任一目标更差则不支配
                if (p.Objectives[i] < q.Objectives[i]) betterInAny = true;
            }
            return betterInAny;
        }

        private static void CalculateCrowdingDistance(List<Individual> front)
        {
            if (front.Count == 0) return;

            foreach (var ind in front)
                ind.CrowdingDistance = 0.0;

            int m = front[0].Objectives.Length;

            for (int obj = 0; obj < m; obj++)
            {
                var sorted = front.OrderBy(ind => ind.Objectives[obj]).ToList();

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

        // ===== 工具 & 日志 =====

        /// <summary>
        /// 从当前种群计算“夏最小 / 冬最大”。注意：我们最小化的是 {夏, -冬}，
        /// 因此原始“冬最大” = -min(f1)。
        /// </summary>
        private static (double bestSummer, double bestWinter) GetGenerationMetrics(List<Individual> pop)
        {
            double bestSummer = pop.Min(ind => ind.Objectives[0]);
            double bestWinter = -pop.Min(ind => ind.Objectives[1]); // 还原为“冬季小时数越大越好”
            return (bestSummer, bestWinter);
        }

        private void ExportFrontToCsv(List<Individual> front, string path)
        {
            using var sw = new StreamWriter(path, false);
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

        private void ExportBestGenes(List<Individual> front, string path)
        {
            var best = front.OrderBy(ind => ind.Objectives.Sum()).First();
            using var sw = new StreamWriter(path, false);
            sw.WriteLine(string.Join(",", best.Genes.Select(v => v.ToString("G17"))));
        }
    }
}
