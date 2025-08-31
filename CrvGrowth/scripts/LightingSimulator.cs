using System;
using System.IO;
using System.Globalization;
using System.Collections.Generic;
using System.Numerics;
using CrvGrowth.Solar; // 需要同项目中的 SolarNoaa.cs

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

        // —— 站点与坐标系（可配置）
        private double _latitudeDeg  = 32.0603;   // 南京
        private double _longitudeDeg = 118.7969;  // 南京
        private double _tzOffsetHours = 8.0;      // UTC+8（不考虑夏令时）

        private Vector3 _up    = new(0, 0, 1);    // Up=+Z
        private Vector3 _north = new(0, 1, 0);    // 北向=+Y（南向外法线=-Y）

        // —— 太阳计算选项
        private bool _useApparentElevation = true; // true=视高度（含折射），false=几何高度
        private double _minElevationDeg = 0.0;     // 小于该高度（含）时跳过（视作无直射）
        
        
        
        
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
            double gridSize)
        {
            if (verticalCurve.Count != extrudedCurve.Count)
                throw new ArgumentException("verticalCurve 和 extrudedCurve 的点数必须相同");

            _verticalCurve = verticalCurve;
            _extrudedCurve = extrudedCurve;
            _date = date;
            _startTime = startTime;
            _endTime = endTime;
            _interval = interval;

            _roomWidth = roomWidth;
            _roomDepth = roomDepth;
            _gridSize = gridSize;

            InitializeGrid();
        }

        private void InitializeGrid()
        {
            _gridCols = (int)Math.Ceiling(_roomWidth / _gridSize);
            _gridRows = (int)Math.Ceiling(_roomDepth / _gridSize);
            _gridCenters = new Vector3[_gridCols, _gridRows];
            _lightHourGrid = new int[_gridCols, _gridRows];

            for (int x = 0; x < _gridCols; x++)
            {
                for (int y = 0; y < _gridRows; y++)
                {
                    float cx = (float)((x + 0.5) * _gridSize);
                    float cy = (float)((y + 0.5) * _gridSize);
                    _gridCenters[x, y] = new Vector3(cx, cy, 0f);
                    _lightHourGrid[x, y] = 0;
                }
            }

            //Console.WriteLine($"初始化网格成功：列数 {_gridCols}, 行数 {_gridRows}, 单元格大小 {_gridSize}");
        }

        public void RunSimulation()
        {
            //Console.WriteLine("开始进行每日多时段光照模拟（含动态太阳角度）...");

            for (var currentTime = _startTime; currentTime <= _endTime; currentTime = currentTime.Add(_interval))
            {
                double hour = currentTime.Hour + currentTime.Minute / 60.0;

                //Vector3 sunDir = Vector3.Normalize(new Vector3(0, 1, -1)); // 或使用 GetSunDirection(hour)
                
                
                
                
                var dtLocal = new DateTime(_date.Year, _date.Month, _date.Day, currentTime.Hour, currentTime.Minute, 0, DateTimeKind.Unspecified);
                
                var angles = SolarNoaa.Compute(
                    dtLocal, _latitudeDeg, _longitudeDeg, _tzOffsetHours,
                    applyRefraction: _useApparentElevation);
                
                // 低于阈值（含地平线）则跳过——无直射
                double el = _useApparentElevation ? angles.ApparentElevationDeg : angles.GeometricElevationDeg;
                if (el <= _minElevationDeg) continue;

                // 世界坐标下“指向太阳”的向量 → 投影用“从太阳指向地面”的方向
                var toSun = SolarNoaa.DirectionToSun(el, angles.AzimuthDeg, _up, _north);
                var sunDir = -toSun;

                // 若几乎切向（Z 分量过小）则跳过，避免数值不稳
                if (Math.Abs(sunDir.Z) < 1e-8) continue;
                
                Console.WriteLine($"→ 当前时刻 {hour}, 太阳角度方向：{sunDir}");
                
                
                

                bool[,] shadowGrid = new bool[_gridCols, _gridRows];

                for (int i = 0; i < _verticalCurve.Count - 1; i++)
                {
                    var v0 = _verticalCurve[i];
                    var v1 = _verticalCurve[i + 1];
                    var b1 = _extrudedCurve[i + 1];
                    var b0 = _extrudedCurve[i];

                    var p0 = ProjectOntoXY(v0, sunDir);
                    var p1 = ProjectOntoXY(v1, sunDir);
                    var p2 = ProjectOntoXY(b1, sunDir);
                    var p3 = ProjectOntoXY(b0, sunDir);

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
                            Vector3 center = _gridCenters[x, y];
                            if (PointInQuad(center, p0, p1, p2, p3))
                            {
                                shadowGrid[x, y] = true;
                            }
                        }
                    }
                }

                for (int x = 0; x < _gridCols; x++)
                {
                    for (int y = 0; y < _gridRows; y++)
                    {
                        if (!shadowGrid[x, y])
                        {
                            _lightHourGrid[x, y] += 1;
                        }
                    }
                }
            }

            //Console.WriteLine("所有时段光照模拟完成，已更新累计光照小时矩阵。");
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

            //Console.WriteLine($"累计光照小时数（两行格式）已保存到：{filePath}");
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

        private double GetSolarAngle(double hour)
        {
            if (hour <= 12.0)
                return 25.0 + 40.0 * ((hour - 8.0) / 4.0);  // 25° → 65°
            else
                return 65.0 - 40.0 * ((hour - 12.0) / 4.0); // 65° → 25°
        }

        private Vector3 GetSunDirection(double hour)
        {
            double thetaDeg = GetSolarAngle(hour);
            double thetaRad = thetaDeg * Math.PI / 180.0;
            float y = (float)Math.Cos(thetaRad);
            float z = (float)-Math.Sin(thetaRad);
            return Vector3.Normalize(new Vector3(0, y, z));
        }
        
        public double GetTotalLightHours()
        {
            double total = 0;
            for (int x = 0; x < _gridCols; x++)
            {
                for (int y = 0; y < _gridRows; y++)
                {
                    total += _lightHourGrid[x, y];
                }
            }
            return total;
        }

        public double GetAverageLightHours()
        {
            double total = GetTotalLightHours();
            return total / (_gridCols * _gridRows);
        }
    }
}
