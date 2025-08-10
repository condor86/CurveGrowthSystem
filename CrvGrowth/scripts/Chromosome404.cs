using GeneticSharp;
using System;
using System.Linq;

namespace CrvGrowth
{
    public class Chromosome404 : ChromosomeBase
    {
        public Chromosome404() : base(404)
        {
            CreateGenes();
        }

        public override Gene GenerateGene(int index)
        {
            // 前 4 个：repeller 影响力，范围 [0.01, 5.00]
            if (index < 4)
            {
                return new Gene(RandomizationProvider.Current.GetDouble(0.01, 5.0));
            }
            // 后 400 个：Z 拉伸量，范围 [0.0, 100.0]
            else
            {
                return new Gene(RandomizationProvider.Current.GetDouble(0.0, 100.0));
            }
        }

        public override IChromosome CreateNew()
        {
            return new Chromosome404();
        }

        public double[] ToArray()
        {
            return GetGenes().Select(g => (double)g.Value).ToArray();
        }
    }
}