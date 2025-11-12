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
    /// LightingSimulator
    /// - 投影：将立面带（verticalCurve + extrudedCurve）按太阳方向投影到 XY；
    /// - 平铺：对投影后的四边形，以 MirrorOffset 为周期在 X/Y 无限平铺；
    /// - 取交：仅枚举与房间 [0,roomWidth)×[0,roomDepth) 有交的副本；
    /// - 栅格：按 GridSize 对房间网格覆盖；未被遮挡的格点每步+1；
    /// - 统计：样本数与小时换算分离（hours = samples × interval.TotalHours）。
    ///
    /// 关键参数：
    ///   RoomWidth/RoomDepth —— 房间尺寸（半开区间），决定网格宽高；
    ///   GridSize            —— 栅格化分辨率（网格涂色单元边长）；
    ///   MirrorOffset        —— 投影后周期平铺的周期（默认 1000）。
    /// </summary>
    public class LightingSimulator
    {
        private readonly List<Vector3> _verticalCurve;
        private readonly List<Vector3> _extrudedCurve;

        private readonly DateOnly _date;
        private readonly TimeOnly _startTime;
        private readonly TimeOnly _endTime;
        private readonly TimeSpan _interval;

        // 房间采用第一象限半开矩形：[0, roomWidth) × [0, roomDepth)
        private readonly double _roomWidth;    // X 尺寸
        private readonly double _roomDepth;    // Y 尺寸

        private readonly double _gridSize;     // 栅格边长（用于涂色/累计）
        private readonly double _mirrorOffset; // 投影后平铺周期（用于复制副本）

        private readonly bool _isClosed;              // true: 末点与首点补段
        private readonly bool _enablePeriodicTiling;  // true: 启用投影后平铺

        // —— 站点与坐标系（NOAA 回退路径）——
        private double _latitudeDeg   = 32.0603;   // 南京
        private double _longitudeDeg  = 118.7969;  // 南京
        private double _tzOffsetHours = 8.0;       // UTC+8（不考虑夏令时）

        private Vector3 _up    = new(0, 0, 1);     // Up=+Z
        private Vector3 _north = new(0, 1, 0);     // 北=+Y（南向外法线=-Y）

        // —— 太阳计算选项 —— 
        private bool   _useApparentElevation = true; // true=视高度，false=几何高度
        private double _minElevationDeg      = 0.0;  // ≤ 阈值则不计直射

        // —— 网格 —— 
        private Vector3[,] _gridCenters;
        private int _gridCols;
        private int _gridRows;

        private int[,] _lightHourGrid;      // 样本计数（每步 +1）
        private bool[,] _shadowGridBuffer;  // 遮挡标记复用缓冲

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
            bool enablePeriodicTiling = true, // 启用“投影后按 MirrorOffset 平铺”
            double mirrorOffset = 1000.0       // 新：平铺周期（默认 1000）
        )
        {
            if (verticalCurve == null || extrudedCurve == null)
                throw new ArgumentNullException("verticalCurve/extrudedCurve 不能为空");
            if (verticalCurve.Count != extrudedCurve.Count)
                throw new ArgumentException("verticalCurve 和 extrudedCurve 的点数必须相同");
            if (roomWidth <= 0 || roomDepth <= 0)
                throw new ArgumentOutOfRangeException("房间尺寸必须为正数");
            if (gridSize <= 0)
                throw new ArgumentOutOfRangeException("gridSize 必须为正数");
            if (mirrorOffset <= 0)
                throw new ArgumentOutOfRangeException("mirrorOffset 必须为正数");

            _verticalCurve = verticalCurve;
            _extrudedCurve = extrudedCurve;

            _date      = date;
            _startTime = startTime;
            _endTime   = endTime;
            _interval  = interval;

            _roomWidth    = roomWidth;
            _roomDepth    = roomDepth;
            _gridSize     = gridSize;
            _mirrorOffset = mirrorOffset;

            _isClosed             = isClosed;
            _enablePeriodicTiling = enablePeriodicTiling;

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
        /// 优化路径：使用已缓存的“指向太阳”的单位向量序列（Z≈sin(高度角)）。
        /// 统一归一化并按高度角阈值过滤。
        /// </summary>
        public void RunWithSunVectors(Vector3[] toSuns)
        {
            if (toSuns == null || toSuns.Length == 0) return;

            double minElSin = Math.Sin(_minElevationDeg * Math.PI / 180.0);

            foreach (var raw in toSuns)
            {
                var norm = raw;
                float len = norm.Length();
                if (len <= 1e-12f) continue;
                norm /= len;

                if (norm.Z <= minElSin) continue; // 低于地平线/阈值不计
                AccumulateForSunVector(norm);
            }
        }

        /// <summary>
        /// 回退路径：按时间步调用 NOAA 计算太阳方位与高度。
        /// </summary>
        public void RunSimulation()
        {
            for (var t = _startTime; t <= _endTime; t = t.Add(_interval))
            {
                var dtLocal = new DateTime(_date.Year, _date.Month, _date.Day,
                                           t.Hour, t.Minute, 0, DateTimeKind.Unspecified);

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
        /// 单步累计：投影 →（可选）平铺取交 → 栅格化 → 未遮挡格点样本+1
        /// </summary>
        private void AccumulateForSunVector(Vector3 toSun)
        {
            var sunDir = -Vector3.Normalize(toSun); // 从太阳指向地面
            if (Math.Abs(sunDir.Z) < 1e-8) return;  // 近切向，数值不稳则跳过

            var shadowGrid = _shadowGridBuffer;
            ClearShadowGrid(shadowGrid);

            int n = _verticalCurve.Count;
            if (n < 2) return;

            // 主段
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

                if (_enablePeriodicTiling)
                    RasterizeQuadTiledIntoRoom(p0, p1, p2, p3, ref shadowGrid);
                else
                    RasterizeQuadToShadowGrid(p0, p1, p2, p3, ref shadowGrid);
            }

            // 闭合补段
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

                if (_enablePeriodicTiling)
                    RasterizeQuadTiledIntoRoom(p0, p1, p2, p3, ref shadowGrid);
                else
                    RasterizeQuadToShadowGrid(p0, p1, p2, p3, ref shadowGrid);
            }

            // 累计：未被遮挡的格点样本 +1
            for (int x = 0; x < _gridCols; x++)
                for (int y = 0; y < _gridRows; y++)
                    if (!shadowGrid[x, y]) _lightHourGrid[x, y] += 1;
        }

        // ====== 平铺与栅格化 ======

        private static Vector3 OffsetXY(in Vector3 p, float dx, float dy)
            => new Vector3(p.X + dx, p.Y + dy, p.Z);

        /// <summary>
        /// 闭式求解：给定 [min,max] 的周期平移 [min + k*s, max + k*s] 与房间 [0,room) 有交的整数 k 范围。
        /// 条件：max + k*s > 0 且 min + k*s < room。
        /// </summary>
        private static (int kmin, int kmax) TileIndexRange(double min, double max, double room, double s)
        {
            if (s <= 0) return (1, 0); // 空集
            const double eps = 1e-12;

            // k > (-max)/s  → kmin = floor((-max)/s + eps) + 1
            int kmin = (int)Math.Floor((-max) / s + eps) + 1;

            // k < (room - min)/s → kmax = floor(((room - min) - eps)/s)
            int kmax = (int)Math.Floor(((room - min) - eps) / s);

            if (kmax < kmin) return (1, 0);
            return (kmin, kmax);
        }

        /// <summary>
        /// 对“投影后四边形”按 MirrorOffset 在 X/Y 方向无限平铺，仅枚举与房间有交的副本并栅格化。
        /// </summary>
        private void RasterizeQuadTiledIntoRoom(
            Vector3 p0, Vector3 p1, Vector3 p2, Vector3 p3, ref bool[,] shadowGrid)
        {
            // 投影四边形 AABB（等价于对四个顶点求并集边界）
            float minX = MathF.Min(MathF.Min(p0.X, p1.X), MathF.Min(p2.X, p3.X));
            float maxX = MathF.Max(MathF.Max(p0.X, p1.X), MathF.Max(p2.X, p3.X));
            float minY = MathF.Min(MathF.Min(p0.Y, p1.Y), MathF.Min(p2.Y, p3.Y));
            float maxY = MathF.Max(MathF.Max(p0.Y, p1.Y), MathF.Max(p2.Y, p3.Y));

            double s = _mirrorOffset; // 使用 MirrorOffset 作为周期

            // 闭式范围（若为空则与房间无交）
            var (kxMin, kxMax) = TileIndexRange(minX, maxX, _roomWidth,  s);
            var (kyMin, kyMax) = TileIndexRange(minY, maxY, _roomDepth, s);
            if (kxMin > kxMax || kyMin > kyMax) return;

            // 枚举所有 (kx, ky)
            for (int kx = kxMin; kx <= kxMax; kx++)
            {
                float dx = (float)(kx * s);
                for (int ky = kyMin; ky <= kyMax; ky++)
                {
                    float dy = (float)(ky * s);

                    var q0 = OffsetXY(p0, dx, dy);
                    var q1 = OffsetXY(p1, dx, dy);
                    var q2 = OffsetXY(p2, dx, dy);
                    var q3 = OffsetXY(p3, dx, dy);

                    RasterizeQuadToShadowGrid(q0, q1, q2, q3, ref shadowGrid);
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

            int minCol = Math.Max(0, (int)Math.Floor((double)minX / _gridSize));
            int minRow = Math.Max(0, (int)Math.Floor((double)minY / _gridSize));

            // 半开区间上界（避免将恰在右/上边界之外的列/行算入）
            int maxCol = Math.Min(_gridCols - 1, (int)Math.Floor(((double)maxX - 1e-7) / _gridSize));
            int maxRow = Math.Min(_gridRows - 1, (int)Math.Floor(((double)maxY - 1e-7) / _gridSize));

            // 完全越界 → 不栅格化
            if (maxCol < 0 || maxRow < 0 || minCol > _gridCols - 1 || minRow > _gridRows - 1)
                return;

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

        // ====== 点测：四边形 → 两三角分解 ======

        private bool PointInQuad(Vector3 p, Vector3 a, Vector3 b, Vector3 c, Vector3 d)
            => PointInTri(p, a, b, c) || PointInTri(p, a, c, d);

        private bool PointInTri(Vector3 p, Vector3 a, Vector3 b, Vector3 c)
        {
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
            // AccumulateForSunVector 已保证 |dir.Z| 足够大
            float t = -p.Z / dir.Z;
            return new Vector3(p.X + t * dir.X, p.Y + t * dir.Y, 0f);
        }

        private void ClearShadowGrid(bool[,] grid) => Array.Clear(grid, 0, grid.Length);

        // ====== 输出与统计 ======

        /// <summary>
        /// 以两行一单元输出：第一行 "{x, y, 0.0}"；第二行样本计数。
        /// </summary>
        public void SaveLightHourGrid(string filePath)
        {
            using StreamWriter writer = new StreamWriter(filePath);
            for (int y = 0; y < _gridRows; y++)
            {
                for (int x = 0; x < _gridCols; x++)
                {
                    var pt = _gridCenters[x, y];
                    int samples = _lightHourGrid[x, y];

                    string coordLine = $"{{{pt.X.ToString(CultureInfo.InvariantCulture)}, " +
                                       $"{pt.Y.ToString(CultureInfo.InvariantCulture)}, 0.0}}";
                    writer.WriteLine(coordLine);
                    writer.WriteLine(samples);
                }
            }
        }

        public double GetTotalLightSamples()
        {
            double total = 0;
            for (int x = 0; x < _gridCols; x++)
                for (int y = 0; y < _gridRows; y++)
                    total += _lightHourGrid[x, y];
            return total;
        }

        public double GetAverageLightSamples()
        {
            double total = GetTotalLightSamples();
            return total / (_gridCols * _gridRows);
        }

        public double GetTotalLightHours()   => GetTotalLightSamples()   * _interval.TotalHours;
        public double GetAverageLightHours() => GetAverageLightSamples() * _interval.TotalHours;
    }
}
