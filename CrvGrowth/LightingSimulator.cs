using System;
using System.IO;
using System.Globalization;
using System.Collections.Generic;
using System.Numerics;
using NumSharp;

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

            Console.WriteLine($"初始化网格成功：列数 {_gridCols}, 行数 {_gridRows}, 单元格大小 {_gridSize}");
        }

        public void RunSimulation()
        {
            Console.WriteLine("开始进行每日多时段光照模拟（含动态太阳角度）...");

            for (var currentTime = _startTime; currentTime <= _endTime; currentTime = currentTime.Add(_interval))
            {
                double hour = currentTime.Hour + currentTime.Minute / 60.0;

                Vector3 sunDir = Vector3.Normalize(new Vector3(0, 1, -1)); // 或使用 GetSunDirection(hour)
                Console.WriteLine($"→ 当前时刻 {hour}, 太阳角度方向：{sunDir}");

                bool[,] shadowGrid = new bool[_gridCols, _gridRows];

                int N = _verticalCurve.Count - 1;
                NDArray quadArray = np.empty(new Shape(N, 4, 3), np.float32);

                for (int i = 0; i < N; i++)
                {
                    quadArray[i, 0] = new float[] { _verticalCurve[i].X,     _verticalCurve[i].Y,     _verticalCurve[i].Z };
                    quadArray[i, 1] = new float[] { _verticalCurve[i + 1].X, _verticalCurve[i + 1].Y, _verticalCurve[i + 1].Z };
                    quadArray[i, 2] = new float[] { _extrudedCurve[i + 1].X, _extrudedCurve[i + 1].Y, _extrudedCurve[i + 1].Z };
                    quadArray[i, 3] = new float[] { _extrudedCurve[i].X,     _extrudedCurve[i].Y,     _extrudedCurve[i].Z };
                }

                NDArray projected = ProjectQuadArrayOntoXY(quadArray, sunDir);  // [N, 4, 3]

                for (int i = 0; i < N; i++)
                {
                    var quad = new Vector3[4];
                    for (int j = 0; j < 4; j++)
                    {
                        var pt = projected[i, j];
                        quad[j] = new Vector3(
                            pt[0].GetSingle(),
                            pt[1].GetSingle(),
                            0f
                        );
                    }

                    Vector3 p0 = quad[0], p1 = quad[1], p2 = quad[2], p3 = quad[3];

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

            Console.WriteLine("所有时段光照模拟完成，已更新累计光照小时矩阵。");
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

            Console.WriteLine($"累计光照小时数（两行格式）已保存到：{filePath}");
        }

        private NDArray ProjectQuadArrayOntoXY(NDArray quads, Vector3 sunDir)
        {
            Vector3 dir = Vector3.Normalize(sunDir);
            float dx = dir.X, dy = dir.Y, dz = dir.Z;

            var reshaped = quads.reshape(-1, 3);  // [N*4, 3]
            var px = reshaped[":", 0];
            var py = reshaped[":", 1];
            var pz = reshaped[":", 2];

            var t = np.negative(pz) / dz;
            var proj_x = px + t * dx;
            var proj_y = py + t * dy;
            var proj_z = np.zeros_like(proj_x);

            var projected = np.stack(new NDArray[] { proj_x, proj_y, proj_z }, axis: 1);  // [N*4, 3]
            return projected.reshape(quads.shape);  // [N, 4, 3]
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
    }
}
