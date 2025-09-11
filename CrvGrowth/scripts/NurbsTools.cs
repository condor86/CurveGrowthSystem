// NurbsTools.cs
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Numerics;

public static class NurbsTools
{
    // 简单的 NURBS 曲线容器（Rhino 默认：w=1，开区间均匀内部结点）
    public sealed class NurbsCurve
    {
        public int Degree;                 // p
        public List<Vector3> Ctrl;         // 控制点 {P_i}
        public double[] Weights;           // 权重 {w_i}（这里全 1）
        public double[] Knots;             // 结点向量 U（非递减）

        public int ControlPointCount => Ctrl.Count; // = n+1
    }

    /// <summary>
    /// 用“Rhino 控制点曲线”的默认方式构造 NURBS：Degree=3、w=1、开区间、内部结点等距
    /// 注：若点数不足，将自动降低次数到 (count-1)
    /// </summary>
    public static NurbsCurve BuildRhinoLikeCurve(IReadOnlyList<Vector3> points, int degree = 3)
    {
        if (points == null || points.Count < 2)
            throw new ArgumentException("需要至少 2 个点来构造曲线");

        int p = Math.Min(degree, points.Count - 1);
        int n = points.Count - 1;                   // 最后一个控制点索引
        int m = n + p + 1;                          // 最后一个结点索引（Knots 长度 = m+1）

        // 结点向量：开区间（两端重复 p+1 次），内部结点等分到 (0,1)
        var U = new double[m + 1];
        for (int i = 0; i <= p; i++) U[i] = 0.0;
        for (int i = m - p; i <= m; i++) U[i] = 1.0;

        int internalCount = n - p;                  // 内部结点个数
        for (int j = 1; j <= internalCount; j++)
            U[p + j] = (double)j / (internalCount + 1);

        // 权重全 1（非有理）
        var W = new double[points.Count];
        for (int i = 0; i < W.Length; i++) W[i] = 1.0;

        return new NurbsCurve
        {
            Degree = p,
            Ctrl = new List<Vector3>(points),
            Weights = W,
            Knots = U
        };
    }

    /// <summary>
    /// 等距弧长采样（返回曲线上均匀弧长间隔的点）
    /// dense 指定密集预采样数量（用于构造 u→弧长 的近似映射）
    /// </summary>
    public static List<Vector3> SampleByArcLength(NurbsCurve c, int count, int dense = -1)
    {
        if (count < 2) throw new ArgumentException("采样点数应 ≥ 2");
        if (dense <= 0) dense = Math.Max(2000, count * 20);

        // 1) 预采样：等参数取样，建立 u 与弧长的离散映射
        var us = new double[dense];
        var pts = new Vector3[dense];
        for (int i = 0; i < dense; i++)
        {
            double u = (double)i / (dense - 1);
            us[i] = u;
            pts[i] = Evaluate(c, u);
        }

        // 2) 累积弧长
        var cum = new double[dense];
        cum[0] = 0.0;
        for (int i = 1; i < dense; i++)
            cum[i] = cum[i - 1] + Vector3.Distance(pts[i - 1], pts[i]);

        double total = cum[dense - 1];

        // 退化情况：曲线长度 ~ 0
        if (total <= 1e-9)
        {
            var flat = new List<Vector3>(count);
            for (int k = 0; k < count; k++) flat.Add(pts[0]);
            return flat;
        }

        // 3) 目标等距长度 → 反查 u，再精确求值
        var outPts = new List<Vector3>(count);
        outPts.Add(pts[0]); // 起点
        for (int k = 1; k < count - 1; k++)
        {
            double target = total * k / (count - 1);
            int idx = LowerBound(cum, target); // 找到 cum[idx] >= target 的最小 idx
            int i0 = Math.Max(0, idx - 1);
            int i1 = idx;

            double seg = cum[i1] - cum[i0];
            double t = seg > 0 ? (target - cum[i0]) / seg : 0.0;

            // 线性插值得到参数 u，再评估曲线（比直接插值点更精确）
            double u = Lerp(us[i0], us[i1], t);
            outPts.Add(Evaluate(c, u));
        }
        outPts.Add(pts[dense - 1]); // 终点
        return outPts;
    }

    /// <summary>
    /// 保存 CSV：每行 x,y,z（InvariantCulture）
    /// </summary>
    public static void SaveCsv(string path, IReadOnlyList<Vector3> points)
    {
        using var sw = new StreamWriter(path);
        var ci = CultureInfo.InvariantCulture;
        for (int i = 0; i < points.Count; i++)
        {
            var p = points[i];
            sw.WriteLine(string.Create(ci, $"{p.X},{p.Y},{p.Z}"));
        }
    }

    // —— 以下为 NURBS 求值所需基础算法（The NURBS Book A2.1/A2.2）——

    /// <summary> 在结点向量 U 中寻找包含 u 的区间跨度下标（span ∈ [p, n]） </summary>
    private static int FindSpan(int n, int p, double u, double[] U)
    {
        // 特判末端
        if (u >= U[n + 1]) return n;
        if (u <= U[p]) return p;

        int low = p, high = n + 1, mid = (low + high) / 2;
        while (u < U[mid] || u >= U[mid + 1])
        {
            if (u < U[mid]) high = mid;
            else low = mid;
            mid = (low + high) / 2;
        }
        return mid;
    }

    /// <summary> 基函数值 N_{i,p}(u)（只返回与 span 相关的 p+1 个非零值） </summary>
    private static void BasisFuns(int span, double u, int p, double[] U, double[] N /*len p+1*/)
    {
        var left = new double[p + 1];
        var right = new double[p + 1];

        N[0] = 1.0;
        for (int j = 1; j <= p; j++)
        {
            left[j] = u - U[span + 1 - j];
            right[j] = U[span + j] - u;
            double saved = 0.0;

            for (int r = 0; r < j; r++)
            {
                double denom = right[r + 1] + left[j - r];
                double temp = denom != 0 ? N[r] / denom : 0.0;
                N[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            N[j] = saved;
        }
    }

    /// <summary> 评估曲线点 C(u)（支持权重） </summary>
    public static Vector3 Evaluate(NurbsCurve c, double u)
    {
        int n = c.ControlPointCount - 1;
        int p = c.Degree;
        var U = c.Knots;
        var W = c.Weights;
        var P = c.Ctrl;

        u = Math.Clamp(u, U[p], U[n + 1]); // 夹到有效参数域

        int span = FindSpan(n, p, u, U);
        var N = new double[p + 1];
        BasisFuns(span, u, p, U, N);

        Vector3 C = Vector3.Zero;
        double wsum = 0.0;

        int start = span - p;
        for (int j = 0; j <= p; j++)
        {
            int i = start + j;
            double w = W[i] * N[j];
            C += (float)w * P[i];
            wsum += w;
        }
        if (wsum == 0.0) return C; // 退化保护
        return (1.0f / (float)wsum) * C;
    }

    // 工具：二分下界 & 线性插值
    private static int LowerBound(double[] arr, double x)
    {
        int lo = 0, hi = arr.Length - 1;
        while (lo < hi)
        {
            int mid = (lo + hi) >> 1;
            if (arr[mid] < x) lo = mid + 1;
            else hi = mid;
        }
        return lo;
    }
    private static double Lerp(double a, double b, double t) => a + (b - a) * t;
}
