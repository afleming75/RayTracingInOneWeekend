
#ifndef RAYTRACER_UTILITIES_H
#define RAYTRACER_UTILITIES_H
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>

const double pi = 3.1415926535897932385;

inline double clamp(double x, double min, double max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}
inline double degrees_to_radians(double degrees) {
    return degrees * pi / 180;
}
//const double infinity = std::numeric_limits<double>::infinity();

inline double ffmin(double a, double b) { return a <= b ? a : b; }
inline double ffmax(double a, double b) { a >= b ? a : b; }

#endif //RAYTRACER_UTILITIES_H
