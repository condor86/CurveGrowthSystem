// File: CrvGrowth/LightingSimulator.cs
using System;
using System.IO;
using System.Globalization;
using System.Collections.Generic;
using System.Numerics;
using CrvGrowth.Solar; // 仅供 RunSimulation() 回退路径使用（NOAA 计算）

namespace CrvGrowth
{
    /// <summary>
    /// 镜像投影模式：
    /// Off      → 仅本块；
    /// LeftRight→ 本块 + 左右（X 方向 ±gridSize）；
    /// Four     → 本块 + 上下左右（X 与 Y 方向各 ±gridSize）。不包含对角（±X±Y）。
    /// </summary>
    public enum MirrorShadowMode
    {
        Off = 0,
        LeftRight = 1,
        Four = 2
    }

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

        private readonly bool _isClosed;                 // 是否按闭合曲线处理（默认 true）
        private readonly MirrorShadowMode _mirrorMode;   // 镜像投影模式

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

        // 复用遮挡网格，避免每步分配
        private bool[,] _shadowGridBuffer;

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
            bool isClosed = true,
            MirrorShadowMode mirrorMode = MirrorShadowMode.LeftRight // 默认保持原行为（左右镜像）
        )
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

            _isClosed   = isClosed;
            _mirrorMode = mirrorMode;

            InitializeGrid();
        }

        private void InitializeGrid()
        {
            _gridCols = (int)Math.Ceiling(_roomWidth / _gridSize);
            _gridRows = (int)Math.Ceiling(_roomDepth / _gridSize);
            _gridCenters      = new Vector3[_gridCols, _gridRows];
            _lightHourGrid    = new int[_gridCols, _gridRows];
            _shadowGridBuffer = new bool[_gridCols, _gridRows];

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

            var shadowGrid = _shadowGridBuffer;
            ClearShadowGrid(shadowGrid);

            int n = _verticalCurve.Count;
            if (n < 2) return;

            // —— 主循环段（开放/闭合同步覆盖 0..n-2）——
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

                RasterizeQuadToShadowGridWithCopies(p0, p1, p2, p3, ref shadowGrid);
            }

            // —— 闭合补段：末尾 → 开头 ——（如 _isClosed）
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

                RasterizeQuadToShadowGridWithCopies(p0, p1, p2, p3, ref shadowGrid);
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

        // ====== 栅格化（含镜像复制开关） ======

        // 平移 XY（保持 Z 不变）
        private static Vector3 OffsetXY(in Vector3 p, float dx, float dy)
        {
            return new Vector3(p.X + dx, p.Y + dy, p.Z);
        }

        /// <summary>
        /// 根据 _mirrorMode 进行栅格化：
        /// - Off：仅原四边形；
        /// - LeftRight：原四边形 + X 方向 ±gridSize；
        /// - Four：原四边形 + X/Y 方向各 ±gridSize（四邻），不包含对角。
        /// </summary>
        private void RasterizeQuadToShadowGridWithCopies(
            Vector3 p0, Vector3 p1, Vector3 p2, Vector3 p3, ref bool[,] shadowGrid)
        {
            // 原位置
            RasterizeQuadToShadowGrid(p0, p1, p2, p3, ref shadowGrid);

            if (_mirrorMode == MirrorShadowMode.Off) return;

            float s = (float)_gridSize;

            if (_mirrorMode == MirrorShadowMode.LeftRight || _mirrorMode == MirrorShadowMode.Four)
            {
                // 左/右（±X）
                RasterizeQuadToShadowGrid(
                    OffsetXY(p0, -s, 0), OffsetXY(p1, -s, 0),
                    OffsetXY(p2, -s, 0), OffsetXY(p3, -s, 0),
                    ref shadowGrid);

                RasterizeQuadToShadowGrid(
                    OffsetXY(p0,  s, 0), OffsetXY(p1,  s, 0),
                    OffsetXY(p2,  s, 0), OffsetXY(p3,  s, 0),
                    ref shadowGrid);
            }

            if (_mirrorMode == MirrorShadowMode.Four)
            {
                // 上/下（±Y）
                RasterizeQuadToShadowGrid(
                    OffsetXY(p0, 0, -s), OffsetXY(p1, 0, -s),
                    OffsetXY(p2, 0, -s), OffsetXY(p3, 0, -s),
                    ref shadowGrid);

                RasterizeQuadToShadowGrid(
                    OffsetXY(p0, 0,  s), OffsetXY(p1, 0,  s),
                    OffsetXY(p2, 0,  s), OffsetXY(p3, 0,  s),
                    ref shadowGrid);
            }

            // 说明：如需对角（±X±Y）可在此追加四个偏移。
        }

        private void RasterizeQuadToShadowGrid(
            Vector3 p0, Vector3 p1, Vector3 p2, Vector3 p3, ref bool[,] shadowGrid)
        {
            float minX = MathF.Min(MathF.Min(p0.X, p1.X), MathF.Min(p2.X, p3.X));
            float maxX = MathF.Max(MathF.Max(p0.X, p1.X), MathF.Max(p2.X, p3.X));
            float minY = MathF.Min(MathF.Min(p0.Y, p1.Y), MathF.Min(p2.Y, p3.Y));
            float maxY = MathF.Max(MathF.Max(p0.Y, p1.Y), MathF.Max(p2.Y, p3.Y));

            int minCol = Math.Max(0, (int)Math.Floor((double)minX / _gridSize));
            int minRow = Math.Max(0, (int)Math.Floor((double)minY / _gridSize));

            // 半开区间上界（避免把恰好在右/上边界之外的列/行算入）
            int maxCol = Math.Min(_gridCols - 1, (int)Math.Floor(((double)maxX - 1e-7) / _gridSize));
            int maxRow = Math.Min(_gridRows - 1, (int)Math.Floor(((double)maxY - 1e-7) / _gridSize));

            if (maxCol < 0 || maxRow < 0 || minCol > _gridCols - 1 || minRow > _gridRows - 1)
                return; // 完全越界

            minCol = Math.Max(0, minCol);
            minRow = Math.Max(0, minRow);
            maxCol = Math.Min(_gridCols - 1, maxCol);
            maxRow = Math.Min(_gridRows - 1, maxRow);

            for (int x = minCol; x <= maxCol; x++)
            {
                for (int y = minRow; y <= maxRow; y++)
                {
                    if (PointInQuad(_gridCenters[x, y], p0, p1, p2, p3))
                        shadowGrid[x, y] = true;
                }
            }
        }

        // ====== 点测：四边形 → 三角分解，更稳健 ======

        private bool PointInQuad(Vector3 p, Vector3 a, Vector3 b, Vector3 c, Vector3 d)
        {
            // 将四边形 a-b-c-d 拆为两个三角形：a-b-c 与 a-c-d
            return PointInTri(p, a, b, c) || PointInTri(p, a, c, d);
        }

        private bool PointInTri(Vector3 p, Vector3 a, Vector3 b, Vector3 c)
        {
            // 基于叉积符号一致性（允许共线视为在内）
            static float Cross(Vector3 u, Vector3 v) => u.X * v.Y - u.Y * v.X;

            var ab = b - a; var ap = p - a;
            var bc = c - b; var bp = p - b;
            var ca = a - c; var cp = p - c;

            float c1 = Cross(ab, ap);
            float c2 = Cross(bc, bp);
            float c3 = Cross(ca, cp);

            const float eps = 1e-6f;
            bool nonNeg = (c1 >= -eps) && (c2 >= -eps) && (c3 >= -eps);
            bool nonPos = (c1 <=  eps) && (c2 <=  eps) && (c3 <=  eps);
            return nonNeg || nonPos;
        }

        // ====== 投影与工具 ======

        private Vector3 ProjectOntoXY(Vector3 p, Vector3 dir)
        {
            // 假设 AccumulateForSunVector 已保证 |dir.Z| 足够大
            float t = -p.Z / dir.Z;
            return new Vector3(p.X + t * dir.X, p.Y + t * dir.Y, 0f);
        }

        private void ClearShadowGrid(bool[,] grid)
        {
            Array.Clear(grid, 0, grid.Length);
        }

        // ====== 输出与统计 ======

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
