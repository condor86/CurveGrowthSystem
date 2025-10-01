// File: scripts/SunCache.cs
using System;
using System.Collections.Generic;
using System.Numerics;
using CrvGrowth.Solar; // 复用你现有的 SolarNoaa.cs

namespace CrvGrowth.Scripts
{
    public static class SunVectors
    {
        private static DateTime At(DateOnly d, TimeOnly t)
            => new DateTime(d.Year, d.Month, d.Day, t.Hour, t.Minute, 0, DateTimeKind.Unspecified);

        /// 生成“指向太阳”的单位向量序列（含起点与终点；不会过滤夜间样本）
        public static Vector3[] Build(
            DateOnly date, TimeOnly start, TimeOnly end, TimeSpan interval,
            double latDeg, double lonDeg, double tzHours,
            Vector3 up, Vector3 north,
            bool useApparentElevation = true)
        {
            var list = new List<Vector3>();
            for (var t = start; t <= end; t = t.Add(interval))
            {
                var dtLocal = At(date, t);
                var ang = SolarNoaa.Compute(dtLocal, latDeg, lonDeg, tzHours, applyRefraction: useApparentElevation);
                var toSun = SolarNoaa.DirectionToSun(
                    useApparentElevation ? ang.ApparentElevationDeg : ang.GeometricElevationDeg,
                    ang.AzimuthDeg, up, north);
                list.Add(toSun);
            }
            return list.ToArray(); // 例如 08,10,12,14,16 → 长度=5
        }
    }
}