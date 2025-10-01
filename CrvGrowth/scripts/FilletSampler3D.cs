using System;
using System.Collections.Generic;
using System.Numerics;

namespace CrvGrowth
{
    public static class FilletUtil
    {
        /// <summary>
        /// 逐拐点圆角化：
        /// - isClosed=false：视为开口折线，保留首尾原点，仅处理中间拐点；
        /// - isClosed=true：把所有点都当作拐点参与圆角，但不做首尾补点的几何闭合（不追加 result[0]）。
        /// 采样方式：按每 90° 的采样密度采样弧。
        /// </summary>
        public static List<Vector3> FilletPolyline(
            IReadOnlyList<Vector3> pts,
            float radius,
            int arcSamplesPer90Deg = 6,
            float angleEpsDeg = 1.0f,
            bool isClosed = false,
            bool clampRadius = true)
        {
            var result = new List<Vector3>();
            if (pts == null || pts.Count < 2) return new List<Vector3>(pts ?? Array.Empty<Vector3>());

            int n = pts.Count;
            bool treatClosed = isClosed && n >= 3;

            if (!treatClosed) result.Add(pts[0]); // closed: 不预先加入首点

            int start = treatClosed ? 0 : 1;
            int end   = treatClosed ? n - 1 : n - 2;

            for (int i = start; i <= end; i++)
            {
                int iPrev = Mod(i - 1, n);
                int iCurr = Mod(i, n);
                int iNext = Mod(i + 1, n);

                var pPrev = pts[iPrev];
                var pCurr = pts[iCurr];
                var pNext = pts[iNext];

                if (!TryFilletCorner(
                        pPrev, pCurr, pNext,
                        radius, Deg2Rad(angleEpsDeg), clampRadius,
                        out var t1, out var t2, out var center, out double arcAngle, out var nrm))
                {
                    // 统一回退策略：失败则保留当前拐点
                    if (result.Count == 0 || !NearlyEqual(result[^1], pCurr))
                        result.Add(pCurr);
                    continue;
                }

                if (result.Count == 0 || !NearlyEqual(result[^1], t1))
                    result.Add(t1);

                foreach (var q in SampleArc(t1, t2, center, nrm, arcAngle, arcSamplesPer90Deg))
                {
                    if (result.Count == 0 || !NearlyEqual(result[^1], q))
                        result.Add(q);
                }
            }

            if (!treatClosed)
            {
                if (result.Count == 0 || !NearlyEqual(result[^1], pts[^1]))
                    result.Add(pts[^1]);
            }
            // closed: 不进行任何首尾闭合的补点操作

            return result;
        }

        /// <summary>
        /// 加强版：每个圆角用固定点数（默认 9，含 T1/T2）替代原中间拐点；
        /// - isClosed=false：视为开口折线，保留首尾原点；
        /// - isClosed=true：把所有点当作拐点参与圆角，但不做首尾补点的几何闭合。
        /// </summary>
        public static List<Vector3> FilletPolylineWithFixedArcPoints(
            IReadOnlyList<Vector3> pts,
            float radius,
            int arcPointCount = 9,
            float angleEpsDeg = 1.0f,
            bool isClosed = false,
            bool clampRadius = true)
        {
            var result = new List<Vector3>();
            if (pts == null || pts.Count < 2) return new List<Vector3>(pts ?? Array.Empty<Vector3>());

            int n = pts.Count;
            bool treatClosed = isClosed && n >= 3;

            if (!treatClosed) result.Add(pts[0]); // closed: 不预先加入首点

            int start = treatClosed ? 0 : 1;
            int end   = treatClosed ? n - 1 : n - 2;

            for (int i = start; i <= end; i++)
            {
                int iPrev = Mod(i - 1, n);
                int iCurr = Mod(i,     n);
                int iNext = Mod(i + 1, n);

                var pPrev = pts[iPrev];
                var pCurr = pts[iCurr];
                var pNext = pts[iNext];

                if (!TryFilletCorner(
                        pPrev, pCurr, pNext,
                        radius, Deg2Rad(angleEpsDeg), clampRadius,
                        out var t1, out var t2, out var center, out double _, out var nrm))
                {
                    // 统一回退策略：失败则保留当前拐点
                    if (result.Count == 0 || !NearlyEqual(result[^1], pCurr))
                        result.Add(pCurr);
                    continue;
                }

                if (result.Count == 0 || !NearlyEqual(result[^1], t1))
                    result.Add(t1);

                foreach (var q in SampleArcFixedCount(t1, t2, center, nrm, arcPointCount))
                {
                    if (result.Count == 0 || !NearlyEqual(result[^1], q))
                        result.Add(q);
                }
            }

            if (!treatClosed)
            {
                if (result.Count == 0 || !NearlyEqual(result[^1], pts[^1]))
                    result.Add(pts[^1]);
            }
            // closed: 不进行任何首尾闭合的补点操作

            return result;
        }

        // —— double 半径重载（便于 const double 直接传入） ——
        public static List<Vector3> FilletPolyline(
            IReadOnlyList<Vector3> pts,
            double radius,
            int arcSamplesPer90Deg = 6,
            float angleEpsDeg = 1.0f,
            bool isClosed = false,
            bool clampRadius = true)
            => FilletPolyline(pts, (float)radius, arcSamplesPer90Deg, angleEpsDeg, isClosed, clampRadius);

        public static List<Vector3> FilletPolylineWithFixedArcPoints(
            IReadOnlyList<Vector3> pts,
            double radius,
            int arcPointCount = 9,
            float angleEpsDeg = 1.0f,
            bool isClosed = false,
            bool clampRadius = true)
            => FilletPolylineWithFixedArcPoints(pts, (float)radius, arcPointCount, angleEpsDeg, isClosed, clampRadius);

        // ===================== 内部实现 =====================

        private static bool TryFilletCorner(
            Vector3 pPrev, Vector3 pCurr, Vector3 pNext,
            float R_in, double angleEpsRad, bool clampRadius,
            out Vector3 T1, out Vector3 T2, out Vector3 O, out double arcAngle, out Vector3 nrm)
        {
            T1 = T2 = O = nrm = default;
            arcAngle = 0.0;

            var a = pCurr - pPrev;
            var b = pNext - pCurr;

            double L1 = a.Length();
            double L2 = b.Length();
            if (L1 < 1e-7 || L2 < 1e-7) return false;

            var v1 = Vector3.Normalize(a);
            var v2 = Vector3.Normalize(b);

            // 局部法向（确定平面与转向）
            var cross = Vector3.Cross(v1, v2);
            double crossLen = cross.Length();
            if (crossLen < 1e-9) return false; // 近共线
            nrm = cross / (float)crossLen;

            // 内角（入射方向 -v1 到 v2）
            double cosTheta = Math.Clamp(Vector3.Dot(-v1, v2), -1f, 1f);
            double theta = Math.Acos(cosTheta); // [0, π]
            if (theta < angleEpsRad) return false; // 近直线

            double R = R_in;
            double half = theta * 0.5;
            double tanHalf = Math.Tan(half);

            // —— 正确公式：沿两边后退距离 d = R / tan(θ/2)（= R * cot(θ/2)）——
            if (Math.Abs(tanHalf) < 1e-12) return false; // 极端退化
            double t = R / tanHalf;

            // 边长约束与半径压缩（若需要）
            double tMax = Math.Min(L1, L2) * 0.999;
            if (t > tMax)
            {
                if (!clampRadius) return false;
                R = tMax * tanHalf;   // 注意：由 t = R / tanHalf 推回 R = t * tanHalf
                t = tMax;
                if (R <= 1e-9) return false;
            }

            // 切点
            T1 = pCurr - v1 * (float)t;
            T2 = pCurr + v2 * (float)t;

            // 圆心：内角平分方向，距离为 R / sin(θ/2)
            double sinHalf = Math.Sin(half);
            if (Math.Abs(sinHalf) < 1e-12) return false;

            var bis = -v1 + v2;
            double bisLen = bis.Length();
            if (bisLen < 1e-9) return false;
            var bdir = bis / (float)bisLen;

            O = pCurr + bdir * (float)(R / sinHalf);

            // 圆心角（取小弧）
            var u1 = Vector3.Normalize(T1 - O);
            var u2 = Vector3.Normalize(T2 - O);
            double ang = SignedAngle(u1, u2, nrm);
            if (ang < 0) ang += 2 * Math.PI;
            if (ang > Math.PI) ang = 2 * Math.PI - ang;
            arcAngle = ang;

            return true;
        }

        /// <summary>固定数量等角度采样（含端点 T1/T2）。</summary>
        private static IEnumerable<Vector3> SampleArcFixedCount(
            Vector3 T1, Vector3 T2, Vector3 O, Vector3 nrm, int count)
        {
            int k = Math.Max(2, count);

            var rVec = T1 - O;
            float r = rVec.Length();
            if (r < 1e-12f)
            {
                yield return T1;
                if (!NearlyEqual(T1, T2)) yield return T2;
                yield break;
            }

            var u1 = rVec / r;
            var v1 = Vector3.Normalize(Vector3.Cross(nrm, u1));

            var u2 = Vector3.Normalize(T2 - O);
            double ang = SignedAngle(u1, u2, nrm);
            if (ang < 0) ang += 2 * Math.PI;
            if (ang > Math.PI) ang = 2 * Math.PI - ang;

            for (int i = 0; i < k; i++)
            {
                double t = (k == 1) ? 0.0 : (double)i / (k - 1);
                double a = ang * t;
                var dir = (float)Math.Cos(a) * u1 + (float)Math.Sin(a) * v1;
                yield return O + dir * r;
            }
        }

        /// <summary>按每 90° 的段数采样（含端点）。</summary>
        private static IEnumerable<Vector3> SampleArc(
            Vector3 T1, Vector3 T2, Vector3 O, Vector3 nrm,
            double arcAngle, int samplesPer90)
        {
            var u1 = Vector3.Normalize(T1 - O);
            var u2 = Vector3.Normalize(T2 - O);
            var u = u1;
            var v = Vector3.Normalize(Vector3.Cross(nrm, u));

            double ang = SignedAngle(u1, u2, nrm);
            if (ang < 0) ang += 2 * Math.PI;
            if (ang > Math.PI) ang = 2 * Math.PI - ang;

            int segs = Math.Max(1, (int)Math.Ceiling(ang / (Math.PI / 2.0)));
            int stepsPerSeg = Math.Max(2, samplesPer90);
            int steps = segs * stepsPerSeg;

            for (int i = 0; i <= steps; i++)
            {
                double t = ang * i / steps;
                var dir = (float)Math.Cos(t) * u + (float)Math.Sin(t) * v;
                var p = O + dir * (T1 - O).Length();
                yield return p;
            }
        }

        private static int Mod(int a, int m) => (a % m + m) % m;
        private static double Deg2Rad(double d) => d * Math.PI / 180.0;

        private static double SignedAngle(Vector3 a, Vector3 b, Vector3 axis)
        {
            var x = Math.Clamp(Vector3.Dot(a, b), -1f, 1f);
            double unsigned = Math.Acos(x);
            var cross = Vector3.Cross(a, b);
            double s = Vector3.Dot(cross, axis);
            return s >= 0 ? unsigned : -unsigned;
        }

        private static bool NearlyEqual(in Vector3 p, in Vector3 q, float eps = 1e-5f)
            => Vector3.DistanceSquared(p, q) <= eps * eps;
    }
}
