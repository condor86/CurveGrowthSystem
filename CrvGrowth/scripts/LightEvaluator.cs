using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace CrvGrowth
{
    public static class LightingEvaluator
    {
        public static (double summer, double winter) Evaluate(
            List<Vector3> baseCurve,   // 约 200 个点
            double[] gene)             // 404 个基因：前4为repeller，后400为extrude量
        {
            int pointCount = baseCurve.Count;
            if (gene.Length < 4 + pointCount)
                throw new ArgumentException("基因数量不足，无法生成完整 extrusion 曲线");

            // 解析基因
            double[] repellerStrengths = gene.Take(4).ToArray();
            double[] zOffsets = gene.Skip(4).Take(pointCount).ToArray();

            // 构建 extruded 曲线（每个点向 -Z 拉伸）
            var extruded = new List<Vector3>(pointCount);
            for (int i = 0; i < pointCount; i++)
            {
                var pt = baseCurve[i];
                var depth = (float)-zOffsets[i];  // 向下拉伸
                extruded.Add(new Vector3(pt.X, pt.Y, pt.Z + depth));
            }

            // 构建 LightingSimulator（夏季）
            var simSummer = new LightingSimulator(
                verticalCurve: baseCurve,
                extrudedCurve: extruded,
                date: new DateOnly(2024, 6, 21),
                startTime: new TimeOnly(8, 0),
                endTime: new TimeOnly(16, 0),
                interval: TimeSpan.FromHours(1),
                roomWidth: 1000,
                roomDepth: 1000,
                gridSize: 10);
            simSummer.RunSimulation();
            double scoreSummer = simSummer.TotalLightHours();

            // 构建 LightingSimulator（冬季）
            var simWinter = new LightingSimulator(
                verticalCurve: baseCurve,
                extrudedCurve: extruded,
                date: new DateOnly(2024, 12, 21),
                startTime: new TimeOnly(8, 0),
                endTime: new TimeOnly(16, 0),
                interval: TimeSpan.FromHours(1),
                roomWidth: 1000,
                roomDepth: 1000,
                gridSize: 10);
            simWinter.RunSimulation();
            double scoreWinter = simWinter.TotalLightHours();

            return (scoreSummer, scoreWinter);  // 注意：主程序中需使用 (scoreSummer, -scoreWinter)
        }
    }
}
