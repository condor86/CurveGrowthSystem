// File: CrvGrowth/Geometry/FilletSampler3D.cs
using System;
using System.Collections.Generic;
using System.Numerics;

namespace CrvGrowth
{
    /// <summary>
    /// 逐顶点 3D 圆角采样（局部平面 fillet）：
    /// 输入：空间折线点集 + 目标半径；输出：仅包含“圆弧采样点 + 端点”的点列。
    /// 直线段不额外采点，由相邻点自然连线。
    /// 采样规则：1/4 圆≈5 点、1/2 圆≈9 点；所有圆弧至少 5 点，更大弧长按比例增密。
    /// </summary>
    public static class FilletSampler3D
    {
        public static List<Vector3> SampleFilletedPolyline(
            List<Vector3> poly,
            double radius = 10.0,
            bool isClosed = false,
            double eps = 1e-7)
        {
            if (poly == null || poly.Count < 2) throw new ArgumentException("poly needs >= 2 points");
            var pts = new List<Vector3>(poly);
            if (isClosed)
            {
                // 已闭合的情况下去掉重复首尾
                if (Vector3.Distance(pts[0], pts[^1]) < (float)eps) pts.RemoveAt(pts.Count - 1);
            }

            int n = pts.Count;
            var outPts = new List<Vector3>(n * 4);

            void Push(Vector3 p)
            {
                if (outPts.Count == 0 || Vector3.Distance(outPts[^1], p) > (float)eps)
                    outPts.Add(p);
            }

            // 开口折线：先推首点
            if (!isClosed) Push(pts[0]);

            int end = isClosed ? n : n - 1;
            for (int i = 0; i < end; i++)
            {
                int iPrev = isClosed ? (i - 1 + n) % n : i - 1;
                int iCurr = i;
                int iNext = isClosed ? (i + 1) % n : i + 1;

                if (!isClosed && (iCurr == 0 || iCurr == n - 1))
                    continue;
                if (iPrev < 0 || iNext >= n) continue;

                Vector3 Pm1 = pts[iPrev];
                Vector3 P   = pts[iCurr];
                Vector3 Pp1 = pts[iNext];

                Vector3 vIn  = Vector3.Normalize(P - Pm1);
                Vector3 vOut = Vector3.Normalize(Pp1 - P);

                float lenIn  = Vector3.Distance(P, Pm1);
                float lenOut = Vector3.Distance(Pp1, P);

                if (lenIn < eps || lenOut < eps)
                {
                    // 极短边：保留原角点，避免不稳定圆角
                    Push(P);
                    continue;
                }

                Vector3 nrm = Vector3.Cross(vIn, vOut);
                if (nrm.Length() < eps)
                {
                    // 近共线：跳过圆角
                    Push(P);
                    continue;
                }
                nrm = Vector3.Normalize(nrm);

                // 在 P 处内角向量
                Vector3 v1 = -vIn; // 指向进入边的反方向
                Vector3 v2 =  vOut;

                double dot = Math.Clamp(Vector3.Dot(v1, v2), -1.0f, 1.0f);
                double alpha = Math.Acos(dot);
                if (alpha < 1e-6 || Math.Abs(alpha - Math.PI) < 1e-6)
                {
                    Push(P);
                    continue;
                }

                // 目标半径与修剪量
                double r = Math.Max(radius, 0.0);
                double t = r * Math.Tan(alpha / 2.0);

                // 半径收缩：保证修剪不超边长
                double tMax = Math.Max(0.0, Math.Min(lenIn, lenOut) - 1e-6);
                if (t > tMax)
                {
                    t = tMax;
                    double tanHalf = Math.Tan(alpha / 2.0);
                    if (tanHalf < 1e-12)
                    {
                        Push(P);
                        continue;
                    }
                    r = t / tanHalf;
                    if (r < 1e-9)
                    {
                        Push(P);
                        continue;
                    }
                }

                // 修剪点
                Vector3 Qin  = P + (Vector3)(v1 * (float)t);
                Vector3 Qout = P + (Vector3)(v2 * (float)t);

                // 圆心：沿角平分线方向
                double sinHalf = Math.Sin(alpha / 2.0);
                if (sinHalf < 1e-12)
                {
                    Push(P);
                    continue;
                }
                double d = r / sinHalf;
                Vector3 bis = Vector3.Normalize(v1 + v2);
                Vector3 C   = P + (Vector3)(bis * (float)d);

                // 局部平面基
                Vector3 e1 = Vector3.Normalize(Qin - C);
                if (e1.Length() < eps)
                {
                    Push(P);
                    continue;
                }
                Vector3 e2 = Vector3.Normalize(Vector3.Cross(nrm, e1));

                static double Atan2InBasis(Vector3 p, Vector3 center, Vector3 e1, Vector3 e2)
                {
                    Vector3 u = p - center;
                    double x = Vector3.Dot(u, e1);
                    double y = Vector3.Dot(u, e2);
                    return Math.Atan2(y, x);
                }

                double angStart = Atan2InBasis(Qin,  C, e1, e2);
                double angEnd   = Atan2InBasis(Qout, C, e1, e2);

                // 转向方向（与法向一致）
                double turnSign = Math.Sign(Vector3.Dot(Vector3.Cross(v1, v2), nrm));
                angEnd = NormalizeAngleToTarget(angStart, angEnd, -(int)turnSign);
    
                double sweep = Math.Abs(angEnd - angStart);
                int sampleCount = ComputeSampleCount(sweep); // 至少 5 点

                // 推入修剪点（避免重复）
                Push(Qin);

                // 中间与末端（包含 Qout）
                for (int k = 1; k < sampleCount; k++)
                {
                    double t01 = (double)k / (sampleCount - 1);
                    double ang = Lerp(angStart, angEnd, t01);
                    Vector3 onArc = C + (Vector3)((float)r * (Vector3)(e1 * (float)Math.Cos(ang) + e2 * (float)Math.Sin(ang)));
                    Push(onArc);
                }
            }

            // 开口折线：推末点
            if (!isClosed) Push(pts[^1]);
            else
            {
                // 闭合：确保首尾重合
                if (Vector3.Distance(outPts[0], outPts[^1]) > (float)eps)
                    outPts.Add(outPts[0]);
            }

            return outPts;
        }

        private static int ComputeSampleCount(double sweep)
        {
            // 目标：θ=π/2 -> 5 点；θ=π -> 9 点；更小弧段也至少 5 点
            double perQuarter = 4.0;
            double baseCount = 1.0 + (sweep / (Math.PI / 2.0)) * perQuarter; // 线性外推
            int count = (int)Math.Round(baseCount, MidpointRounding.AwayFromZero);
            if (count < 5) count = 5;
            return count;
        }

        private static double NormalizeAngleToTarget(double a0, double a1, int dirSign)
        {
            double twoPi = Math.PI * 2.0;
            double delta = a1 - a0;
            delta = delta - Math.Round(delta / twoPi) * twoPi;

            if (dirSign >= 0)
            {
                if (delta < 0) delta += twoPi;  // 正向最短
            }
            else
            {
                if (delta > 0) delta -= twoPi;  // 反向最短
            }
            return a0 + delta;
        }

        private static double Lerp(double a, double b, double t) => a + (b - a) * t;
    }
}
