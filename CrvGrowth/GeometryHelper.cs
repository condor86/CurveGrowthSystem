using System;
using System.Numerics;

namespace CrvGrowth
{
    public static class GeometryHelper
    {
        /// <summary>
        /// 将点 P 沿方向 D 投影到以点 O 为基准、法向量为 N 的平面上。
        /// </summary>
        public static Vector3 ProjectPointOntoPlane(Vector3 P, Vector3 D, Vector3 O, Vector3 N)
        {
            float dDotN = Vector3.Dot(D, N);
            if (Math.Abs(dDotN) < 1e-8)
                throw new InvalidOperationException("投影方向与平面平行，无法投影。");

            Vector3 PO = O - P;
            float t = Vector3.Dot(PO, N) / dDotN;

            return P + t * D;
        }
    }
}