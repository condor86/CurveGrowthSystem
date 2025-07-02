using NetTopologySuite.Geometries;
using System;

namespace CrvGrowth
{
    public readonly struct Point3D(double x, double y, double z = 0.0)
    {
        public double X { get; } = x;
        public double Y { get; } = y;
        public double Z { get; } = z;

        public double DistanceTo(Point3D other)
        {
            double dx = X - other.X;
            double dy = Y - other.Y;
            double dz = Z - other.Z;
            return Math.Sqrt(dx * dx + dy * dy + dz * dz);
        }

        // 运算符重载
        public static Vector3D operator -(Point3D a, Point3D b)
        {
            return new Vector3D(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
        }

        public static Point3D operator +(Point3D p, Vector3D v)
        {
            return new Point3D(p.X + v.X, p.Y + v.Y, p.Z + v.Z);
        }

        public static Point3D operator +(Point3D p, Point3D v)
        {
            return new Point3D(p.X + v.X, p.Y + v.Y, p.Z + v.Z);
        }

        public static Point3D operator *(double t, Point3D p)
        {
            return new Point3D(t * p.X, t * p.Y, t * p.Z);
        }

        public static Point3D Lerp(Point3D a, Point3D b, double t)
        {
            return new Point3D(
                a.X + (b.X - a.X) * t,
                a.Y + (b.Y - a.Y) * t,
                a.Z + (b.Z - a.Z) * t
            );
        }

        public Point3D RotateXYToXZ()
        {
            return new Point3D(X, 0.0, Y); // Y 上升为 Z，Y 被清零
        }
        
        public Envelope ToEnvelope(double buffer = 0.001)
        {
            // 仅考虑 XY 平面用于二维索引
            return new Envelope(X - buffer, X + buffer, Y - buffer, Y + buffer);
        }

        public override string ToString()
        {
            return $"({X:F3}, {Y:F3}, {Z:F3})";
        }
    }
}