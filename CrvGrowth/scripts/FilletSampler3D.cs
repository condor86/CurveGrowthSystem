// File: CrvGrowth/Scripts/FilletSampler3D.cs
using System;
using System.Collections.Generic;
using System.Numerics;

namespace CrvGrowth.Scripts
{
    /// <summary>
    /// 逐顶点 3D 圆角采样（局部平面 fillet）。
    /// - 半径 r 作为固定值（仅当放不下时被动缩小）。
    /// - 张角越大，圆弧扫角越小（≈ π − θ），视觉不鼓包。
    /// - 仅采圆弧点（含圆弧起点/终点=切点与中间点），直线段由相邻点自然直连。
    /// - 采样：1/4 圆≈5 点，1/2 圆≈9 点，所有圆弧至少 5 点。
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
                // 去掉重复首尾
                if (Vector3.Distance(pts[0], pts[^1]) < (float)eps)
                    pts.RemoveAt(pts.Count - 1);
            }

            int n = pts.Count;
            var outPts = new List<Vector3>(n * 4);

            void Push(Vector3 p)
            {
                if (outPts.Count == 0 || Vector3.Distance(outPts[^1], p) > (float)eps)
                    outPts.Add(p);
            }

            // 开口折线：先压首点
            if (!isClosed) Push(pts[0]);

            int end = isClosed ? n : n - 1;
            for (int i = 0; i < end; i++)
            {
                int iPrev = isClosed ? (i - 1 + n) % n : i - 1;
                int iCurr = i;
                int iNext = isClosed ? (i + 1) % n : i + 1;

                if (!isClosed && (iCurr == 0 || iCurr == n - 1)) continue;
                if (iPrev < 0 || iNext >= n) continue;

                Vector3 Pm1 = pts[iPrev];
                Vector3 P   = pts[iCurr];
                Vector3 Pp1 = pts[iNext];

                Vector3 a = P - Pm1; float lenIn  = a.Length();
                Vector3 b = Pp1 - P; float lenOut = b.Length();
                if (lenIn < eps || lenOut < eps) { Push(P); continue; }

                a /= lenIn; // 进入边方向（Pm1→P）
                b /= lenOut; // 离开边方向（P→Pp1）

                // 法向（用于确定局部平面与旋转方向）
                Vector3 nrm = Vector3.Cross(a, b);
                float nrmLen = nrm.Length();
                if (nrmLen < eps)
                {
                    // 共线或近共线：跳过圆角
                    Push(P);
                    continue;
                }
                nrm /= nrmLen;

                // 以 v1 = -a（从 P 指向“进入边”内部），v2 = b（从 P 沿“离开边”）
                Vector3 v1 = -a;
                Vector3 v2 =  b;

                // 内角 α 与转角 θ
                double dot = Math.Clamp(Vector3.Dot(v1, v2), -1.0f, 1.0f);
                double alpha = Math.Acos(dot);      // 内角 ∈ [0, π]
                double theta = Math.PI - alpha;     // 转角（标准 fillet 用它）
                // θ→0：近直线；θ→π：回折，tan(θ/2)→∞，不做圆角
                if (theta < 1e-6 || theta > Math.PI - 1e-6) { Push(P); continue; }

                // 固定半径策略：r 取输入；若放不下再被动缩小
                double tanHalf = Math.Tan(theta / 2.0);
                if (tanHalf < 1e-12) { Push(P); continue; }

                double r = Math.Max(radius, 0.0);
                double t = r * tanHalf; // 每侧修剪量

                double tLimit = Math.Max(0.0, Math.Min(lenIn, lenOut) - 1e-6);
                if (t > tLimit)
                {
                    t = tLimit;
                    r = t / tanHalf;                 // 被动缩小半径（只在放不下时）
                    if (r < 1e-9) { Push(P); continue; }
                }

                // 切点（弧端）—— 与两侧直线相切
                Vector3 Qin  = P + v1 * (float)t;    // 进入边方向上退 t
                Vector3 Qout = P + v2 * (float)t;    // 离开边方向上退 t

                // 圆心：沿角平分线（v1+v2）方向，距离 d = r / sin(θ/2)
                Vector3 bis = v1 + v2;
                float bisLen = bis.Length();
                if (bisLen < eps) { Push(P); continue; } // 数值防护
                bis /= bisLen;

                double sinHalf = Math.Sin(theta / 2.0);
                if (sinHalf < 1e-12) { Push(P); continue; }
                double d = r / sinHalf;

                Vector3 C = P + bis * (float)d;

                // 局部正交基（右手系）：e1 指向 Qin 径向，e2 = e1 × nrm
                Vector3 e1 = Qin - C; float e1Len = e1.Length();
                if (e1Len < eps) { Push(P); continue; }
                e1 /= e1Len;

                Vector3 e2 = Vector3.Cross(e1, nrm);
                float e2Len = e2.Length();
                if (e2Len < eps) { Push(P); continue; }
                e2 /= e2Len;

                static double Atan2InBasis(Vector3 p, Vector3 center, Vector3 e1, Vector3 e2)
                {
                    Vector3 u = p - center;
                    double x = Vector3.Dot(u, e1);
                    double y = Vector3.Dot(u, e2);
                    return Math.Atan2(y, x);
                }

                double angStart = Atan2InBasis(Qin,  C, e1, e2);
                double angEnd   = Atan2InBasis(Qout, C, e1, e2);

                // 期望转向：从 v1 → v2 围绕 nrm 的符号
                int dirSign = Math.Sign(Vector3.Dot(Vector3.Cross(v1, v2), nrm)); // +1 逆时针，-1 顺时针（右手系）

                angEnd = NormalizeAngleToTarget(angStart, angEnd, -dirSign);

                // 计算实际扫角并尽量贴近 θ（必要时尝试±2π与反向修正）
                double sweep = Math.Abs(angEnd - angStart);
                double twoPi = Math.PI * 2.0;
                double diff  = Math.Abs(sweep - theta);

                // 尝试加减 2π 贴合 θ
                double alt1 = Math.Abs(Math.Abs(angEnd + dirSign * twoPi - angStart) - theta);
                if (alt1 + 1e-6 < diff)
                {
                    angEnd += dirSign * twoPi;
                    sweep = Math.Abs(angEnd - angStart);
                    diff  = Math.Abs(sweep - theta);
                }
                // 若仍相差较大，尝试翻转方向（极少见的数值分支）
                if (diff > 1e-3)
                {
                    angEnd = NormalizeAngleToTarget(angStart, angEnd, -dirSign);
                    sweep  = Math.Abs(angEnd - angStart);
                    diff   = Math.Abs(sweep - theta);
                }

                // —— 采样输出 —— //
                // 只采圆弧：先压切点 Qin，再压中间点，最后压“准确的 Qout”（保证端点切向）
                Push(Qin);

                int sampleCount = ComputeSampleCount(sweep); // ≥5
                if (sampleCount > 2)
                {
                    for (int k = 1; k < sampleCount - 1; k++)
                    {
                        double t01 = (double)k / (sampleCount - 1);
                        double ang = Lerp(angStart, angEnd, t01);
                        Vector3 onArc = C + (Vector3)((float)r * (e1 * (float)Math.Cos(ang) + e2 * (float)Math.Sin(ang)));
                        Push(onArc);
                    }
                }

                // 末端切点：用 Qout（而不是参数化计算的终点），确保数值上与离开边切向
                Push(Qout);
            }

            // 开口折线：压尾点；闭合则闭合首尾
            if (!isClosed) Push(pts[^1]);
            else
            {
                if (Vector3.Distance(outPts[0], outPts[^1]) > (float)eps)
                    outPts.Add(outPts[0]);
            }

            return outPts;
        }

        private static int ComputeSampleCount(double sweep)
        {
            // θ=π/2 -> 5 点；θ=π -> 9 点；更小弧也至少 5 点
            double perQuarter = 4.0;
            double baseCount = 1.0 + (sweep / (Math.PI / 2.0)) * perQuarter;
            int count = (int)Math.Round(baseCount, MidpointRounding.AwayFromZero);
            if (count < 5) count = 5;
            return count;
        }

        private static double NormalizeAngleToTarget(double a0, double a1, int dirSign)
        {
            double twoPi = Math.PI * 2.0;
            // 规范到 (-π, π]
            double delta = a1 - a0;
            delta = delta - Math.Round(delta / twoPi) * twoPi;

            if (dirSign >= 0)
            {
                if (delta < 0) delta += twoPi;   // 取正向最短
            }
            else
            {
                if (delta > 0) delta -= twoPi;   // 取反向最短
            }
            return a0 + delta;
        }

        private static double Lerp(double a, double b, double t) => a + (b - a) * t;
    }
}
