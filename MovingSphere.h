#ifndef RAYTRACER_MOVINGSPHERE_H
#define RAYTRACER_MOVINGSPHERE_H
#include <memory>
#include "Vec3.h"
#include "Ray.h"
#include "Material.h"
#include <memory>

using std::shared_ptr;

class MovingSphere: public Hittable {
public:
    MovingSphere() {}
    MovingSphere( Vec3 cen0, Vec3 cen1, double t0, double t1, double r, shared_ptr<Material> mat_ptr) :
            center0(cen0), center1(cen1), time0(t0), time1(t1), radius(r), mat_ptr(mat_ptr) {};

    virtual bool hit(const Ray& r, double tmin, double tmax, HitRecord& rec) const;
    virtual bool bounding_box(double t0, double t1, Aabb& output_box) const;
    Vec3 center(double time) const;

    Vec3 center0, center1;
    double time0, time1;
    double radius;
    shared_ptr<Material> mat_ptr;
};

bool MovingSphere::bounding_box(double t0, double t1, Aabb& output_box) const {
    Aabb box0(center(t0) - Vec3(radius, radius, radius),
              center(t0) + Vec3(radius, radius, radius));
    Aabb box1(center(t1) - Vec3(radius, radius, radius),
              center(t1) + Vec3(radius, radius, radius));
    output_box = surrounding_box(box0, box1);
    return true;
}

Vec3 MovingSphere::center(double time) const {
    return center0 + ((time - time0) / (time1 - time0))*(center1 - center0);
}

bool MovingSphere::hit(const Ray& r, double tmin, double tmax, HitRecord& rec) const {
    Vec3 oc = r.get_origin() - center(r.get_time());
    auto a = r.get_direction().squared_length();
    auto half_b = dot(oc, r.get_direction()); //b = 2half_b
    auto c = oc.squared_length() - radius*radius;
    auto discriminant = half_b*half_b - a*c;
    if (discriminant > 0) {
        auto temp = (-half_b - sqrt(discriminant)) / a;
        if (temp < tmax && temp > tmin) {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            Vec3 outward_normal = (rec.p - center(r.get_time())) / radius;
            rec.normal = outward_normal;
            rec.mat_ptr = mat_ptr;
            return true;
        }
        temp = (-half_b + sqrt(discriminant)) / a;
        if (temp < tmax && temp > tmin) {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            Vec3 outward_normal = (rec.p - center(r.get_time())) / radius;
            rec.normal = outward_normal;
            rec.mat_ptr = mat_ptr;
            return true;
        }
    }
    return false;
}
#endif //RAYTRACER_MOVINGSPHERE_H
