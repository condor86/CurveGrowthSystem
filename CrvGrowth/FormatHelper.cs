// FormattingHelper.cs
using System.Globalization;
using NumSharp;

namespace CrvGrowth
{
    public static class FormattingHelper
    {
        /// <summary>
        /// 格式化 NDArray 点为 {x, y, z} 形式，保留高精度小数。
        /// </summary>
        public static string FormatNDArrayXYZ(NDArray pt)
        {
            return $"{{{((double)pt[0]).ToString(CultureInfo.InvariantCulture)}, " +
                   $"{((double)pt[1]).ToString(CultureInfo.InvariantCulture)}, " +
                   $"{((double)pt[2]).ToString(CultureInfo.InvariantCulture)}}}";
        }

        /// <summary>
        /// 格式化 NDArray 点为 {x, y, 0.0}（专用于光照网格结果输出）。
        /// </summary>
        public static string FormatNDArrayXY0(NDArray pt)
        {
            return $"{{{((double)pt[0]).ToString(CultureInfo.InvariantCulture)}, " +
                   $"{((double)pt[1]).ToString(CultureInfo.InvariantCulture)}, 0.0}}";
        }
    }
}