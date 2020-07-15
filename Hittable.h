
#ifndef RAYTRACER_HITTABLE_H
#define RAYTRACER_HITTABLE_H
#include "Ray.h"
#include "Aabb.h"

#include <memory>
using std::shared_ptr;

class Material;

struct HitRecord {
    double t;
    double u;
    double v;
    Vec3 p;
    Vec3 normal;
    shared_ptr<Material> mat_ptr;
    bool front_face;

    inline void set_face_normal(const Ray& r, const Vec3& outward_normal) {
        front_face = dot(r.get_direction(), outward_normal) < 0;
        normal = front_face ? outward_normal :-outward_normal;
    }
};

class Hittable {
public:
    virtual bool hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const = 0;
    virtual bool bounding_box(double t0, double t1, Aabb& output_box) const = 0;
};
#endif //RAYTRACER_HITTABLE_H
