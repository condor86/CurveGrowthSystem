using System;
using System.Collections.Generic;
using System.Linq;

namespace NSGAII
{
    public class Individual
    {
        public double[] Genes;         // 染色体
        public double[] Objectives;    // 多目标值
        public int Rank;               // 非支配等级
        public double CrowdingDistance;

        public Individual(int geneLength)
        {
            Genes = new double[geneLength];
            Objectives = new double[2];  // 默认两个目标
        }

        public Individual Clone()
        {
            var copy = new Individual(Genes.Length);
            Genes.CopyTo(copy.Genes, 0);
            Objectives.CopyTo(copy.Objectives, 0);
            copy.Rank = Rank;
            copy.CrowdingDistance = CrowdingDistance;
            return copy;
        }
    }

    public class NSGAII
    {
        public int PopulationSize = 100;
        public int Generations = 50;
        public int GeneLength = 10;
        public double CrossoverRate = 0.9;
        public double MutationRate = 0.05;
        public Random Rand = new Random();

        public List<Individual> Run()
        {
            var population = InitializePopulation();

            for (int gen = 0; gen < Generations; gen++)
            {
                EvaluatePopulation(population);
                var offspring = Reproduce(population);
                EvaluatePopulation(offspring);

                var combined = population.Concat(offspring).ToList();
                var sorted = FastNonDominatedSort(combined);

                population = SelectNewPopulation(sorted);
                Console.WriteLine($"Generation {gen + 1}: Best Front Size = {sorted[0].Count}");
            }

            return population;
        }

        private List<Individual> InitializePopulation()
        {
            var list = new List<Individual>();
            for (int i = 0; i < PopulationSize; i++)
            {
                var ind = new Individual(GeneLength);
                for (int j = 0; j < GeneLength; j++)
                    ind.Genes[j] = Rand.NextDouble(); // 默认范围 [0, 1]
                list.Add(ind);
            }
            return list;
        }

        private void EvaluatePopulation(List<Individual> pop)
        {
            foreach (var ind in pop)
            {
                // 示例：两个目标分别是 gene 平均值 与 方差（可自定义）
                var avg = ind.Genes.Average();
                var var = ind.Genes.Select(g => (g - avg) * (g - avg)).Average();

                ind.Objectives[0] = avg;
                ind.Objectives[1] = var;
            }
        }

        private List<Individual> Reproduce(List<Individual> pop)
        {
            var offspring = new List<Individual>();
            while (offspring.Count < PopulationSize)
            {
                var p1 = TournamentSelect(pop);
                var p2 = TournamentSelect(pop);

                Individual c1, c2;
                Crossover(p1, p2, out c1, out c2);
                Mutate(c1);
                Mutate(c2);

                offspring.Add(c1);
                if (offspring.Count < PopulationSize) offspring.Add(c2);
            }
            return offspring;
        }

        private Individual TournamentSelect(List<Individual> pop)
        {
            var a = pop[Rand.Next(pop.Count)];
            var b = pop[Rand.Next(pop.Count)];
            if (a.Rank < b.Rank) return a;
            if (a.Rank > b.Rank) return b;
            return a.CrowdingDistance > b.CrowdingDistance ? a : b;
        }

        private void Crossover(Individual p1, Individual p2, out Individual c1, out Individual c2)
        {
            c1 = p1.Clone();
            c2 = p2.Clone();
            if (Rand.NextDouble() < CrossoverRate)
            {
                int point = Rand.Next(1, GeneLength);
                for (int i = 0; i < point; i++)
                {
                    c1.Genes[i] = p2.Genes[i];
                    c2.Genes[i] = p1.Genes[i];
                }
            }
        }

        private void Mutate(Individual ind)
        {
            for (int i = 0; i < GeneLength; i++)
            {
                if (Rand.NextDouble() < MutationRate)
                    ind.Genes[i] += Rand.NextDouble() * 0.1 - 0.05;  // 微扰
            }
        }

        private List<List<Individual>> FastNonDominatedSort(List<Individual> pop)
        {
            var fronts = new List<List<Individual>>();
            var S = new Dictionary<Individual, List<Individual>>();
            var n = new Dictionary<Individual, int>();
            var rank = new Dictionary<Individual, int>();

            var F1 = new List<Individual>();

            foreach (var p in pop)
            {
                S[p] = new List<Individual>();
                n[p] = 0;

                foreach (var q in pop)
                {
                    if (Dominates(p, q))
                        S[p].Add(q);
                    else if (Dominates(q, p))
                        n[p]++;
                }

                if (n[p] == 0)
                {
                    rank[p] = 0;
                    p.Rank = 0;
                    F1.Add(p);
                }
            }

            fronts.Add(F1);
            int i = 0;
            while (fronts[i].Count > 0)
            {
                var nextFront = new List<Individual>();
                foreach (var p in fronts[i])
                {
                    foreach (var q in S[p])
                    {
                        n[q]--;
                        if (n[q] == 0)
                        {
                            q.Rank = i + 1;
                            rank[q] = i + 1;
                            nextFront.Add(q);
                        }
                    }
                }
                i++;
                fronts.Add(nextFront);
            }

            return fronts;
        }

        private bool Dominates(Individual p, Individual q)
        {
            bool betterInAny = false;
            for (int i = 0; i < p.Objectives.Length; i++)
            {
                if (p.Objectives[i] > q.Objectives[i]) return false;
                if (p.Objectives[i] < q.Objectives[i]) betterInAny = true;
            }
            return betterInAny;
        }

        private List<Individual> SelectNewPopulation(List<List<Individual>> fronts)
        {
            var newPop = new List<Individual>();
            int i = 0;

            while (newPop.Count + fronts[i].Count <= PopulationSize)
            {
                CalculateCrowdingDistance(fronts[i]);
                newPop.AddRange(fronts[i]);
                i++;
            }

            var remaining = PopulationSize - newPop.Count;
            if (remaining > 0)
            {
                CalculateCrowdingDistance(fronts[i]);
                var sorted = fronts[i].OrderByDescending(ind => ind.CrowdingDistance).ToList();
                newPop.AddRange(sorted.Take(remaining));
            }

            return newPop;
        }

        private void CalculateCrowdingDistance(List<Individual> front)
        {
            int nObj = front[0].Objectives.Length;
            foreach (var ind in front)
                ind.CrowdingDistance = 0;

            for (int m = 0; m < nObj; m++)
            {
                var sorted = front.OrderBy(ind => ind.Objectives[m]).ToList();
                sorted[0].CrowdingDistance = double.PositiveInfinity;
                sorted[^1].CrowdingDistance = double.PositiveInfinity;

                double min = sorted[0].Objectives[m];
                double max = sorted[^1].Objectives[m];
                double range = max - min;

                if (range == 0) continue;

                for (int i = 1; i < sorted.Count - 1; i++)
                {
                    sorted[i].CrowdingDistance += (sorted[i + 1].Objectives[m] - sorted[i - 1].Objectives[m]) / range;
                }
            }
        }
    }
}
