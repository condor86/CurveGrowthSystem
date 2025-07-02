using System;

namespace CrvGrowth
{
    public readonly struct Vector3D(double x, double y, double z = 0.0)
    {
        public double X { get; } = x;
        public double Y { get; } = y;
        public double Z { get; } = z;

        public double Length => Math.Sqrt(X * X + Y * Y + Z * Z);

        public static Vector3D Zero => new Vector3D(0, 0, 0);

        public Vector3D Normalized
        {
            get
            {
                double len = Length;
                return len < 1e-8 ? new Vector3D(0, 0, 0) : new Vector3D(X / len, Y / len, Z / len);
            }
        }

        // 运算符重载
        public static Vector3D operator +(Vector3D a, Vector3D b)
        {
            return new Vector3D(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
        }

        public static Vector3D operator -(Vector3D a, Vector3D b)
        {
            return new Vector3D(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
        }

        public static Vector3D operator *(Vector3D v, double scalar)
        {
            return new Vector3D(v.X * scalar, v.Y * scalar, v.Z * scalar);
        }

        public static Vector3D operator *(double scalar, Vector3D v)
        {
            return new Vector3D(v.X * scalar, v.Y * scalar, v.Z * scalar);
        }

        public static Vector3D operator /(Vector3D v, double scalar)
        {
            return new Vector3D(v.X / scalar, v.Y / scalar, v.Z / scalar);
        }

        public override string ToString()
        {
            return $"<{X:F3}, {Y:F3}, {Z:F3}>";
        }
    }
}