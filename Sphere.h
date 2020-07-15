
#ifndef RAYTRACER_SPHERE_H
#define RAYTRACER_SPHERE_H
#include "Hittable.h"
#include "Material.h"
#include "Aabb.h"

class Sphere: public Hittable {
public:
    Sphere() : 	radius(0), mat_ptr(0) {}
    Sphere(Vec3 cen, double r, shared_ptr<Material> mat_ptr) : center(cen), radius(r), mat_ptr(mat_ptr) {}
    virtual bool hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const;
    virtual bool bounding_box(double t0, double t1, Aabb& output_box) const;

    Vec3 center;
    double radius;
    shared_ptr<Material> mat_ptr;
};


bool Sphere::bounding_box(double t0, double t1, Aabb& output_box) const {
    output_box = Aabb(center - Vec3(radius, radius, radius),
                      center + Vec3(radius, radius, radius));
    return true;
}

bool Sphere::hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const {
    Vec3 oc = r.get_origin() - center;
    double a = dot(r.get_direction(), r.get_direction());
    double b = dot(oc, r.get_direction());
    double c = dot(oc, oc) - radius*radius;
    double discriminant = b*b - a*c;
    if (discriminant > 0) {
        double temp = (-b - sqrt(discriminant)) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            Vec3 outward_normal = (rec.p - center) / radius;
            //
            rec.set_face_normal(r, outward_normal);
            rec.mat_ptr = mat_ptr;
            return true;
        }
        temp = (-b + sqrt(discriminant)) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            Vec3 outward_normal = (rec.p - center) / radius;
            //rec.mat_ptr = mat_ptr;
            rec.set_face_normal(r, outward_normal);
            rec.mat_ptr = mat_ptr;
            return true;
        }
    }
    return false;
}

#endif //RAYTRACER_SPHERE_H
