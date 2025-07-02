namespace CrvGrowth
{
    public static class GeometryHelper
    {
        /// <summary>
        /// 将点 P 沿方向 D 投影到以点 O 为基准、法向量为 N 的平面上。
        /// </summary>
        public static Point3D ProjectPointOntoPlane(Point3D P, Vector3D D, Point3D O, Vector3D N)
        {
            double dDotN = D.X * N.X + D.Y * N.Y + D.Z * N.Z;
            if (Math.Abs(dDotN) < 1e-8)
                throw new InvalidOperationException("投影方向与平面平行，无法投影。");

            Vector3D PO = new Vector3D(O.X - P.X, O.Y - P.Y, O.Z - P.Z);
            double t = (PO.X * N.X + PO.Y * N.Y + PO.Z * N.Z) / dDotN;

            return new Point3D(
                P.X + t * D.X,
                P.Y + t * D.Y,
                P.Z + t * D.Z
            );
        }
    }
}