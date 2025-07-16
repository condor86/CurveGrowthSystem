using System.Globalization;
using System.Collections.Generic;
using System.IO;
using NumSharp;

namespace CrvGrowth
{
    public static class IOHelper
    {
        /// <summary>
        /// 从形如 {x, y, z} 的文本文件中加载所有点，返回 NDArray，shape = [N, 3]
        /// </summary>
        public static NDArray LoadPointsAsNDArray(string filePath)
        {
            var rows = new List<double[]>();

            foreach (var line in File.ReadAllLines(filePath))
            {
                var trimmed = line.Trim('{', '}', ' ', '\t');
                var parts = trimmed.Split(',');

                if (parts.Length >= 2 &&
                    double.TryParse(parts[0], NumberStyles.Float, CultureInfo.InvariantCulture, out double x) &&
                    double.TryParse(parts[1], NumberStyles.Float, CultureInfo.InvariantCulture, out double y))
                {
                    double z = 0.0;
                    if (parts.Length >= 3 &&
                        double.TryParse(parts[2], NumberStyles.Float, CultureInfo.InvariantCulture, out double parsedZ))
                        z = parsedZ;

                    rows.Add(new double[] { x, y, z });
                }
            }

            return np.array(rows.ToArray());  // shape: [N, 3]
        }

        /// <summary>
        /// 从单列数字文件中读取 repeller 的 factor，返回 NDArray，shape = [N]
        /// </summary>
        public static NDArray LoadFactorsAsNDArray(string filePath)
        {
            var values = new List<double>();

            foreach (var line in File.ReadAllLines(filePath))
            {
                if (double.TryParse(line.Trim(), NumberStyles.Float, CultureInfo.InvariantCulture, out double val))
                    values.Add(val);
            }

            return np.array(values.ToArray());  // shape: [N]
        }

        /// <summary>
        /// 将 [N, 3] 的 NDArray 保存为文本文件，每行格式为 {x, y, z}
        /// </summary>
        public static void SaveNDArrayAsPointFile(string filePath, NDArray nd)
        {
            using StreamWriter writer = new StreamWriter(filePath);
            int N = nd.shape[0];

            for (int i = 0; i < N; i++)
            {
                string line = "{" +
                    ((double)nd[i, 0]).ToString(CultureInfo.InvariantCulture) + ", " +
                    ((double)nd[i, 1]).ToString(CultureInfo.InvariantCulture) + ", " +
                    ((double)nd[i, 2]).ToString(CultureInfo.InvariantCulture) + "}";

                if (i < N - 1)
                    writer.WriteLine(line);
                else
                    writer.Write(line);  // 最后一行不加换行
            }
        }
    }
}
