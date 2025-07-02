using System.Globalization;

namespace CrvGrowth
{
    public static class IOHelper
    {
        public static List<Point3D> LoadPointsFromFile(string filePath)
        {
            var points = new List<Point3D>();
            foreach (var line in File.ReadAllLines(filePath))
            {
                var trimmed = line.Trim('{', '}', ' ');
                var parts = trimmed.Split(',');
                if (parts.Length >= 2 &&
                    double.TryParse(parts[0], NumberStyles.Float, CultureInfo.InvariantCulture, out double x) &&
                    double.TryParse(parts[1], NumberStyles.Float, CultureInfo.InvariantCulture, out double y))
                {
                    double z = 0.0;
                    if (parts.Length >= 3 &&
                        double.TryParse(parts[2], NumberStyles.Float, CultureInfo.InvariantCulture, out double parsedZ))
                        z = parsedZ;

                    points.Add(new Point3D(x, y, z));
                }
            }
            return points;
        }

        public static List<double> LoadFactorsFromFile(string filePath)
        {
            var factors = new List<double>();
            foreach (var line in File.ReadAllLines(filePath))
            {
                if (double.TryParse(line, NumberStyles.Float, CultureInfo.InvariantCulture, out double f))
                {
                    factors.Add(f);
                }
            }
            return factors;
        }

        public static void SavePointsToFile(string filePath, List<Point3D> points)
        {
            using StreamWriter writer = new StreamWriter(filePath);
            for (int i = 0; i < points.Count; i++)
            {
                var pt = points[i];
                string line = $"{{{pt.X.ToString(CultureInfo.InvariantCulture)}, " +
                              $"{pt.Y.ToString(CultureInfo.InvariantCulture)}, " +
                              $"{pt.Z.ToString(CultureInfo.InvariantCulture)}}}";
                if (i < points.Count - 1)
                    writer.WriteLine(line);
                else
                    writer.Write(line); // 最后一行不加换行
            }
        }
    }
}