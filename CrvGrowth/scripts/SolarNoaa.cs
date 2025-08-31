// SolarNoaa.cs
// NOAA 低阶太阳位置算法（与 Ladybug/LBT 一致的思路）：
// - 方位角定义：自正北 0° 起顺时针（90°=东，180°=南，270°=西）
// - Elevation：提供几何高度与“视高度”（含简式折射修正）
// - 可将 (Elevation, Azimuth) 映射为世界坐标单位向量（指向太阳）
// 仅依赖 System 与 System.Numerics

using System;
using System.Numerics;

namespace CrvGrowth.Solar
{
    public readonly struct SolarAngles
    {
        public readonly double GeometricElevationDeg; // 几何高度角
        public readonly double ApparentElevationDeg;  // 视高度角（含折射修正）
        public readonly double AzimuthDeg;            // 方位角：自正北顺时针 [0,360)
        public readonly double DeclinationDeg;        // 赤纬
        public readonly double HourAngleDeg;          // 时角（真太阳时相对正午）
        public readonly double EquationOfTimeMin;     // 时间方程（分钟）
        public readonly DateTime SolarNoonLocal;      // 当日真太阳时正午（本地时间）

        public SolarAngles(
            double geomEl, double appEl, double az,
            double decl, double hra, double eot, DateTime noonLocal)
        {
            GeometricElevationDeg = geomEl;
            ApparentElevationDeg  = appEl;
            AzimuthDeg            = az;
            DeclinationDeg        = decl;
            HourAngleDeg          = hra;
            EquationOfTimeMin     = eot;
            SolarNoonLocal        = noonLocal;
        }
    }

    public static class SolarNoaa
    {
        /// <summary>
        /// 计算给定“本地时间/纬度/经度/时区”下的太阳角度（NOAA 简式）。
        /// azimuth 从正北 0° 顺时针。applyRefraction=true 时返回“视高度”。
        /// </summary>
        public static SolarAngles Compute(
            DateTime localTime,
            double latitudeDeg, double longitudeDeg, double tzOffsetHours,
            bool applyRefraction = true)
        {
            // —— 本地 → UTC
            if (localTime.Kind == DateTimeKind.Unspecified)
                localTime = DateTime.SpecifyKind(localTime, DateTimeKind.Unspecified);
            var dto = new DateTimeOffset(localTime, TimeSpan.FromHours(tzOffsetHours));
            var utc = dto.ToUniversalTime().DateTime;

            int n = utc.DayOfYear;
            double hourUtc = utc.Hour + utc.Minute / 60.0 + utc.Second / 3600.0;

            // 分数年 γ（弧度）
            double gamma = 2.0 * Math.PI / 365.0 * (n - 1 + (hourUtc - 12.0) / 24.0);

            // 时间方程 EOT（分钟）
            double eot = 229.18 * (0.000075
                         + 0.001868 * Math.Cos(gamma)
                         - 0.032077 * Math.Sin(gamma)
                         - 0.014615 * Math.Cos(2 * gamma)
                         - 0.040849 * Math.Sin(2 * gamma));

            // 太阳赤纬 δ（弧度）
            double delta = 0.006918
                         - 0.399912 * Math.Cos(gamma)
                         + 0.070257 * Math.Sin(gamma)
                         - 0.006758 * Math.Cos(2 * gamma)
                         + 0.000907 * Math.Sin(2 * gamma)
                         - 0.002697 * Math.Cos(3 * gamma)
                         + 0.00148  * Math.Sin(3 * gamma);

            // 真太阳时（分钟）
            double timeOffsetMin = eot + 4.0 * longitudeDeg - 60.0 * tzOffsetHours;
            double trueSolarTimeMin = (localTime.Hour * 60.0 + localTime.Minute + localTime.Second / 60.0 + timeOffsetMin) % 1440.0;

            // 时角 HRA（度；-180~+180）
            double hraDeg = trueSolarTimeMin / 4.0 - 180.0;
            double hraRad = Deg2Rad(hraDeg);

            double phi = Deg2Rad(latitudeDeg);

            // 天顶角 / 几何高度角
            double cosTheta = Math.Sin(phi) * Math.Sin(delta) + Math.Cos(phi) * Math.Cos(delta) * Math.Cos(hraRad);
            cosTheta = Math.Clamp(cosTheta, -1.0, 1.0);
            double theta = Math.Acos(cosTheta);
            double elGeom = 90.0 - Rad2Deg(theta);

            // 折射修正（NOAA 简式）
            double elApp = elGeom;
            if (applyRefraction && elGeom > -0.575)
            {
                double R = 1.02 / Math.Tan(Deg2Rad(elGeom + 10.3 / (elGeom + 5.11))) / 60.0; // arcmin -> deg
                elApp += R;
            }

            // 方位角（自北顺时针）
            double azRad = Math.Atan2(Math.Sin(hraRad),
                             Math.Cos(hraRad) * Math.Sin(phi) - Math.Tan(delta) * Math.Cos(phi));
            double azDeg = (Rad2Deg(azRad) + 180.0) % 360.0;

            // 当日真太阳时正午（本地时间估计）
            var solarNoonMin = SolarNoonLocalMinutes(localTime.Date, longitudeDeg, tzOffsetHours);
            var noonLocal = localTime.Date.AddMinutes(solarNoonMin);

            return new SolarAngles(
                elGeom, elApp, azDeg,
                Rad2Deg(delta), hraDeg, eot, noonLocal);
        }

        /// <summary>
        /// 将 (Elevation, Azimuth) 转为世界坐标单位向量（指向太阳）。
        /// 需提供 Up 与 North（会自动正交化，East=North×Up）。
        /// </summary>
        public static Vector3 DirectionToSun(double elevationDeg, double azimuthDeg, Vector3 up, Vector3 north)
        {
            // 正交化
            up = Vector3.Normalize(up);
            var northProj = north - Vector3.Dot(north, up) * up;
            north = Vector3.Normalize(northProj);
            var east = Vector3.Normalize(Vector3.Cross(north, up)); // 右手系：Y×Z = X

            double el = Deg2Rad(elevationDeg);
            double az = Deg2Rad(azimuthDeg);

            var ground = (float)Math.Cos(az) * north + (float)Math.Sin(az) * east;
            var dir = (float)Math.Sin(el) * up + (float)Math.Cos(el) * ground;
            return Vector3.Normalize(dir); // 指向太阳（从场点朝向太阳）
        }

        // —— 内部辅助 ——
        private static double Deg2Rad(double d) => d * Math.PI / 180.0;
        private static double Rad2Deg(double r) => r * 180.0 / Math.PI;

        private static double SolarNoonLocalMinutes(DateTime dayLocal, double lonDeg, double tzHours)
        {
            // 简式迭代两次足够
            double EstMinutes(DateTime local)
            {
                var dto = new DateTimeOffset(local, TimeSpan.FromHours(tzHours));
                var utc = dto.ToUniversalTime().DateTime;
                int n = utc.DayOfYear;
                double hourUtc = utc.Hour + utc.Minute / 60.0 + utc.Second / 3600.0;
                double gamma = 2.0 * Math.PI / 365.0 * (n - 1 + (hourUtc - 12.0) / 24.0);
                double eot = 229.18 * (0.000075
                             + 0.001868 * Math.Cos(gamma)
                             - 0.032077 * Math.Sin(gamma)
                             - 0.014615 * Math.Cos(2 * gamma)
                             - 0.040849 * Math.Sin(2 * gamma));
                return 720.0 - 4.0 * lonDeg - eot + 60.0 * tzHours;
            }
            double m = 12 * 60;
            m = EstMinutes(dayLocal.AddMinutes(m));
            m = EstMinutes(dayLocal.AddMinutes(m));
            return m;
        }
    }
}
