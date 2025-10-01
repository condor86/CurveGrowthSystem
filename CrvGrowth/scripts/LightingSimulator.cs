using System;
using System.IO;
using System.Globalization;
using System.Collections.Generic;
using System.Numerics;
using CrvGrowth.Solar; // 仅供 RunSimulation() 回退路径使用（NOAA 计算）

namespace CrvGrowth
{
    public class LightingSimulator
    {
        private readonly List<Vector3> _verticalCurve;
        private readonly List<Vector3> _extrudedCurve;

        private readonly DateOnly _date;
        private readonly TimeOnly _startTime;
        private readonly TimeOnly _endTime;
        private readonly TimeSpan _interval;

        private readonly double _roomWidth;
        private readonly double _roomDepth;
        private readonly double _gridSize;

        private readonly bool _isClosed;  // 新增：是否按闭合曲线处理（默认 true）

        // —— 站点与坐标系（用于 NOAA 回退路径）——
        private double _latitudeDeg   = 32.0603;   // 南京
        private double _longitudeDeg  = 118.7969;  // 南京
        private double _tzOffsetHours = 8.0;       // UTC+8（不考虑夏令时）

        private Vector3 _up    = new(0, 0, 1);     // Up=+Z
        private Vector3 _north = new(0, 1, 0);     // 北向=+Y（南向外法线=-Y）

        // —— 太阳计算选项（与 NOAA 回退路径相关）——
        private bool   _useApparentElevation = true; // true=视高度（含折射），false=几何高度
        private double _minElevationDeg      = 0.0;  // ≤该高度视作无直射（含地平线）

        private Vector3[,] _gridCenters;
        private int _gridCols;
        private int _gridRows;

        private int[,] _lightHourGrid;

        public LightingSimulator(
            List<Vector3> verticalCurve,
            List<Vector3> extrudedCurve,
            DateOnly date,
            TimeOnly startTime,
            TimeOnly endTime,
            TimeSpan interval,
            double roomWidth,
            double roomDepth,
            double gridSize,
            bool isClosed = true) // 新增参数：是否闭合（默认 true）
        {
            if (verticalCurve.Count != extrudedCurve.Count)
                throw new ArgumentException("verticalCurve 和 extrudedCurve 的点数必须相同");

            _verticalCurve = verticalCurve;
            _extrudedCurve = extrudedCurve;

            _date      = date;
            _startTime = startTime;
            _endTime   = endTime;
            _interval  = interval;

            _roomWidth = roomWidth;
            _roomDepth = roomDepth;
            _gridSize  = gridSize;

            _isClosed = isClosed;

            InitializeGrid();
        }

        private void InitializeGrid()
        {
            _gridCols = (int)Math.Ceiling(_roomWidth / _gridSize);
            _gridRows = (int)Math.Ceiling(_roomDepth / _gridSize);
            _gridCenters   = new Vector3[_gridCols, _gridRows];
            _lightHourGrid = new int[_gridCols, _gridRows];

            for (int x = 0; x < _gridCols; x++)
            {
                for (int y = 0; y < _gridRows; y++)
                {
                    float cx = (float)((x + 0.5) * _gridSize);
                    float cy = (float)((y + 0.5) * _gridSize);
                    _gridCenters[x, y]   = new Vector3(cx, cy, 0f);
                    _lightHourGrid[x, y] = 0;
                }
            }
        }

        /// <summary>
        /// 优化阶段推荐路径：使用预先计算好的“指向太阳”的单位向量序列，避免在热路径中做 NOAA 计算
        /// </summary>
        public void RunWithSunVectors(Vector3[] toSuns)
        {
            if (toSuns == null || toSuns.Length == 0) return;

            // 用 toSun.Z（=sin(elevation)）与阈值比较，近似地平线过滤
            double minElSin = Math.Sin(_minElevationDeg * Math.PI / 180.0);

            foreach (var toSun in toSuns)
            {
                if (toSun.Z <= minElSin) continue; // 低于地平线/阈值不计
                AccumulateForSunVector(toSun);
            }
        }

        /// <summary>
        /// 回退路径：保持原本的 NOAA 按时刻计算（导出或对比时可用）
        /// </summary>
        public void RunSimulation()
        {
            for (var currentTime = _startTime; currentTime <= _endTime; currentTime = currentTime.Add(_interval))
            {
                var dtLocal = new DateTime(_date.Year, _date.Month, _date.Day,
                                           currentTime.Hour, currentTime.Minute, 0,
                                           DateTimeKind.Unspecified);

                var angles = SolarNoaa.Compute(
                    dtLocal, _latitudeDeg, _longitudeDeg, _tzOffsetHours,
                    applyRefraction: _useApparentElevation);

                double el = _useApparentElevation ? angles.ApparentElevationDeg : angles.GeometricElevationDeg;
                if (el <= _minElevationDeg) continue;

                var toSun = SolarNoaa.DirectionToSun(el, angles.AzimuthDeg, _up, _north);
                AccumulateForSunVector(toSun);
            }
        }

        /// <summary>
        /// 单步累计：给定“指向太阳”的单位向量，完成投影、栅格覆盖与累计
        /// </summary>
        private void AccumulateForSunVector(Vector3 toSun)
        {
            var sunDir = -Vector3.Normalize(toSun); // 从太阳指向地面
            if (Math.Abs(sunDir.Z) < 1e-8) return;  // 近切向，数值不稳则跳过

            bool[,] shadowGrid = new bool[_gridCols, _gridRows];

            int n = _verticalCurve.Count;
            if (n < 2) return;

            // —— 开放/闭合的“主循环段” ——（开放时：0..n-2；闭合时同样先覆盖 0..n-2）
            for (int i = 0; i < n - 1; i++)
            {
                var v0 = _verticalCurve[i];
                var v1 = _verticalCurve[i + 1];
                var b1 = _extrudedCurve[i + 1];
                var b0 = _extrudedCurve[i];

                var p0 = ProjectOntoXY(v0, sunDir);
                var p1 = ProjectOntoXY(v1, sunDir);
                var p2 = ProjectOntoXY(b1, sunDir);
                var p3 = ProjectOntoXY(b0, sunDir);

                RasterizeQuadToShadowGrid(p0, p1, p2, p3, ref shadowGrid);
            }

            // —— 闭合补段：末尾 → 开头 ——（新增）
            if (_isClosed && n >= 2)
            {
                int last = n - 1;
                var v0 = _verticalCurve[last];
                var v1 = _verticalCurve[0];
                var b1 = _extrudedCurve[0];
                var b0 = _extrudedCurve[last];

                var p0 = ProjectOntoXY(v0, sunDir);
                var p1 = ProjectOntoXY(v1, sunDir);
                var p2 = ProjectOntoXY(b1, sunDir);
                var p3 = ProjectOntoXY(b0, sunDir);

                RasterizeQuadToShadowGrid(p0, p1, p2, p3, ref shadowGrid);
            }

            // —— 将未被遮挡的格点累计“光照次数” ——（每个样本步 +1）
            for (int x = 0; x < _gridCols; x++)
            {
                for (int y = 0; y < _gridRows; y++)
                {
                    if (!shadowGrid[x, y])
                        _lightHourGrid[x, y] += 1;
                }
            }
        }

        private void RasterizeQuadToShadowGrid(
            Vector3 p0, Vector3 p1, Vector3 p2, Vector3 p3, ref bool[,] shadowGrid)
        {
            float minX = MathF.Min(MathF.Min(p0.X, p1.X), MathF.Min(p2.X, p3.X));
            float maxX = MathF.Max(MathF.Max(p0.X, p1.X), MathF.Max(p2.X, p3.X));
            float minY = MathF.Min(MathF.Min(p0.Y, p1.Y), MathF.Min(p2.Y, p3.Y));
            float maxY = MathF.Max(MathF.Max(p0.Y, p1.Y), MathF.Max(p2.Y, p3.Y));

            int minCol = Math.Max(0, (int)Math.Floor(minX / _gridSize));
            int maxCol = Math.Min(_gridCols - 1, (int)Math.Ceiling(maxX / _gridSize));
            int minRow = Math.Max(0, (int)Math.Floor(minY / _gridSize));
            int maxRow = Math.Min(_gridRows - 1, (int)Math.Ceiling(maxY / _gridSize));

            for (int x = minCol; x <= maxCol; x++)
            {
                for (int y = minRow; y <= maxRow; y++)
                {
                    if (PointInQuad(_gridCenters[x, y], p0, p1, p2, p3))
                        shadowGrid[x, y] = true;
                }
            }
        }

        public void SaveLightHourGrid(string filePath)
        {
            using StreamWriter writer = new StreamWriter(filePath);
            for (int y = 0; y < _gridRows; y++)
            {
                for (int x = 0; x < _gridCols; x++)
                {
                    var pt = _gridCenters[x, y];
                    int hours = _lightHourGrid[x, y];

                    string coordLine = $"{{{pt.X.ToString(CultureInfo.InvariantCulture)}, " +
                                       $"{pt.Y.ToString(CultureInfo.InvariantCulture)}, 0.0}}";
                    writer.WriteLine(coordLine);
                    writer.WriteLine(hours);
                }
            }
        }

        private Vector3 ProjectOntoXY(Vector3 p, Vector3 dir)
        {
            float t = -p.Z / dir.Z;
            return new Vector3(p.X + t * dir.X, p.Y + t * dir.Y, 0f);
        }

        private bool PointInQuad(Vector3 p, Vector3 a, Vector3 b, Vector3 c, Vector3 d)
        {
            bool SameSide(Vector3 p1, Vector3 p2, Vector3 a1, Vector3 a2)
            {
                float cp1 = (a2.X - a1.X) * (p1.Y - a1.Y) - (a2.Y - a1.Y) * (p1.X - a1.X);
                float cp2 = (a2.X - a1.X) * (p2.Y - a1.Y) - (a2.Y - a1.Y) * (p2.X - a1.X);
                return cp1 * cp2 >= 0;
            }

            return SameSide(p, c, a, b) &&
                   SameSide(p, d, b, c) &&
                   SameSide(p, a, c, d) &&
                   SameSide(p, b, d, a);
        }

        public double GetTotalLightHours()
        {
            double total = 0;
            for (int x = 0; x < _gridCols; x++)
                for (int y = 0; y < _gridRows; y++)
                    total += _lightHourGrid[x, y];
            return total;
        }

        public double GetAverageLightHours()
        {
            double total = GetTotalLightHours();
            return total / (_gridCols * _gridRows);
        }
    }
}
