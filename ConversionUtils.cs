using System;
using static System.Math;

namespace ConsoleGpsConv {
    class ConversionUtils {
    
        //Converto un punto GPS geodetico WGS-84 (latitudine, longitudine, altitudine) in una coordinata ECEF (x, y, z)
        public static void GeodeticToEcef(double lat, double lon, double h, out double x, out double y, out double z) {

            //converto in radianti i gradi sessagesimali/decimali
            var lambda = DegreesToRadians(lat);
            var phi = DegreesToRadians(lon);

            var s = Sin(lambda);
            var N = a / Sqrt(1 - e_sq * s * s);

            var sin_lambda = Sin(lambda);
            var cos_lambda = Cos(lambda);
            var cos_phi = Cos(phi);
            var sin_phi = Sin(phi);

            x = (h + N) * cos_lambda * cos_phi;
            y = (h + N) * cos_lambda * sin_phi;
            z = (h + (1 - e_sq) * N) * sin_lambda;
        }
        
        
        private static double DegreesToRadians(double degrees) {
            return PI / 180.0 * degrees;
        }
        
        
    }
  
}
