#ifndef RAYTRACER_AABB_H
#define RAYTRACER_AABB_H
#include <algorithm>
#include "Vec3.h"
#include "Ray.h"

class Aabb {
public:
    //virtual ~Aabb() {}
    Aabb() {}
    Aabb(const Vec3& a, const Vec3& b) { min = a; max = b; }

    Vec3 get_min() const {return min; }
    Vec3 get_max() const {return max; }

    bool hit(const Ray& r, double tmin, double tmax) const;
    Aabb surrounding_box(Aabb box0, Aabb box1) const;

    Vec3 min, max;
};

Aabb surrounding_box(Aabb box0, Aabb box1) {
    Vec3 small(fmin(box0.get_min().x(), box1.get_min().x()),
               fmin(box0.get_min().y(), box1.get_min().y()),
               fmin(box0.get_min().z(), box1.get_min().z()));
    Vec3 big(fmax(box0.get_max().x(), box1.get_max().x()),
             fmax(box0.get_max().y(), box1.get_max().y()),
             fmax(box0.get_max().z(), box1.get_max().z()));
    return Aabb(small, big);
}

bool Aabb::hit(const Ray& r, double tmin, double tmax) const {
    for(int a = 0; a < 3; a++) {
        auto invD = 1.0f / r.get_direction()[a];
        auto t0 = (get_min()[a] - r.get_origin()[a]) * invD;
        auto t1 = (get_min()[a] - r.get_origin()[a]) * invD;
        if (invD < 0.0f)
            std::swap(t0, t1);
        tmin = t0 > tmin ? t0: tmin;
        tmax = t1 < tmax ? t1: tmax;
        if (tmax <= tmin)
            return false;
    }
    return true;
}
#endif //RAYTRACER_AABB_H
