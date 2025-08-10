using GeneticSharp;
using System;
using System.Collections.Generic;
using System.Numerics;

namespace CrvGrowth
{
    public class FitnessEvaluator : IFitness
    {
        private readonly List<Vector3> _baseCurve;

        public FitnessEvaluator(List<Vector3> baseCurve)
        {
            _baseCurve = baseCurve;
        }

        public double Evaluate(IChromosome chromosome)
        {
            var ch = (Chromosome404)chromosome;
            var genes = ch.ToArray();

            var (summer, winter) = LightingEvaluator.Evaluate(_baseCurve, genes);

            // 注意：NSGA-II 支持多目标，但 IFitness.Evaluate 返回的是单一值
            // 我们将在主程序中使用 IObjectiveEvaluator<double[]> 接口注入两个目标
            // 所以此函数仅供 Debug 时使用
            return -summer + winter; // 不用于真正评价，仅 placeholder
        }
    }
}