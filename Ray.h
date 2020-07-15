#ifndef RAYTRACER_RAY_H
#define RAYTRACER_RAY_H
#include "Vec3.h"

class Ray {
public:
    Vec3 origin;
    Vec3 direction;
    double time;
    Ray() { origin = Vec3(0, 0, 0); direction = Vec3(0, 0, 0); time = 0.0; }
    Ray(const Vec3& orig, const Vec3& direct, double tm = 0.0) : origin(orig), direction(direct), time(tm) {}
    Vec3 get_origin() const { return origin; }
    Vec3 get_direction() const { return direction; }
    double get_time() const { return time; } //used for motion blur, to get average for camera and object movement with one ray calculation
    Vec3 point_at_parameter(double t) const { return origin + t*direction; }
};


#endif //RAYTRACER_RAY_H
