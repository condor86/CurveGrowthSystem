using System.Globalization;

namespace CrvGrowth
{
    public class LightingSimulator
    {
        private readonly List<Point3D> _verticalCurve;
        private readonly List<Point3D> _extrudedCurve;

        private readonly DateOnly _date;
        private readonly TimeOnly _startTime;
        private readonly TimeOnly _endTime;
        private readonly TimeSpan _interval;

        private readonly double _roomWidth;
        private readonly double _roomDepth;
        private readonly double _gridSize;

        private Point3D[,] _gridCenters;
        private int _gridCols;
        private int _gridRows;

        private int[,] _lightHourGrid;

        public LightingSimulator(
            List<Point3D> verticalCurve,
            List<Point3D> extrudedCurve,
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
            _gridCenters = new Point3D[_gridCols, _gridRows];
            _lightHourGrid = new int[_gridCols, _gridRows];

            for (int x = 0; x < _gridCols; x++)
            {
                for (int y = 0; y < _gridRows; y++)
                {
                    double cx = (x + 0.5) * _gridSize;
                    double cy = (y + 0.5) * _gridSize;
                    _gridCenters[x, y] = new Point3D(cx, cy, 0.0);
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

//                Vector3D sunDir = GetSunDirection(hour);
                Vector3D sunDir = new Vector3D(0, 1, -1).Normalized;

                Console.WriteLine($"→ 当前时刻 {hour}, 太阳角度方向：{sunDir}");

                // 每个时间点初始化 shadowGrid
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

                    double minX = Math.Min(Math.Min(p0.X, p1.X), Math.Min(p2.X, p3.X));
                    double maxX = Math.Max(Math.Max(p0.X, p1.X), Math.Max(p2.X, p3.X));
                    double minY = Math.Min(Math.Min(p0.Y, p1.Y), Math.Min(p2.Y, p3.Y));
                    double maxY = Math.Max(Math.Max(p0.Y, p1.Y), Math.Max(p2.Y, p3.Y));

                    int minCol = Math.Max(0, (int)Math.Floor(minX / _gridSize));
                    int maxCol = Math.Min(_gridCols - 1, (int)Math.Ceiling(maxX / _gridSize));
                    int minRow = Math.Max(0, (int)Math.Floor(minY / _gridSize));
                    int maxRow = Math.Min(_gridRows - 1, (int)Math.Ceiling(maxY / _gridSize));

                    for (int x = minCol; x <= maxCol; x++)
                    {
                        for (int y = minRow; y <= maxRow; y++)
                        {
                            Point3D center = _gridCenters[x, y];
                            if (PointInQuad(center, p0, p1, p2, p3))
                            {
                                shadowGrid[x, y] = true;
                            }
                        }
                    }
                }

                // 将当前 shadowGrid 的结果累计到 light hour 统计矩阵中
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

        private Point3D ProjectOntoXY(Point3D p, Vector3D dir)
        {
            double t = -p.Z / dir.Z;
            return new Point3D(p.X + t * dir.X, p.Y + t * dir.Y, 0.0);
        }

        private bool PointInQuad(Point3D p, Point3D a, Point3D b, Point3D c, Point3D d)
        {
            bool SameSide(Point3D p1, Point3D p2, Point3D a1, Point3D a2)
            {
                double cp1 = (a2.X - a1.X) * (p1.Y - a1.Y) - (a2.Y - a1.Y) * (p1.X - a1.X);
                double cp2 = (a2.X - a1.X) * (p2.Y - a1.Y) - (a2.Y - a1.Y) * (p2.X - a1.X);
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

        private Vector3D GetSunDirection(double hour)
        {
            double thetaDeg = GetSolarAngle(hour);
            double thetaRad = thetaDeg * Math.PI / 180.0;
            double y = Math.Cos(thetaRad);
            double z = -Math.Sin(thetaRad);
            return new Vector3D(0, y, z).Normalized;
        }

    }
}
