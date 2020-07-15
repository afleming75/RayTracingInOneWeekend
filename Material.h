
#ifndef RAYTRACER_MATERIAL_H
#define RAYTRACER_MATERIAL_H
#include "Hittable.h"
#include "Ray.h"
#include "Utilities.h"
#include "Texture.h"

double schlick(double cosine, double ref_idx) {
    double r0 = (1 - ref_idx) / (1 + ref_idx);
    r0 = r0*r0;
    return r0 + (1 - r0)*pow((1 - cosine), 5);
}

/*bool refract(const Vec3& v, const Vec3& n, double ni_over_nt, Vec3& refracted) {
    Vec3 uv = unit_vector(v);
    double dt = dot(uv, n);
    double discriminant = 1.0 - ni_over_nt*(1-dt*dt);
    if (discriminant > 0) {
        refracted = ni_over_nt*(uv - n*dt) - n*sqrt(discriminant);
        return true;
    }
    else
        return false;
}*/

Vec3 refract(const Vec3& uv, const Vec3& n, double etai_over_etat) {
    auto cos_theta = dot(-uv, n);
    Vec3 r_out_parallel =  etai_over_etat * (uv + cos_theta*n);
    Vec3 r_out_perp = -sqrt(1.0 - r_out_parallel.squared_length()) * n;
    return r_out_parallel + r_out_perp;
}

Vec3 reflect(const Vec3& v, const Vec3& n) {
    return v - 2*dot(v, n)*n;
}

/*Vec3 random_in_unit_sphere() {
    Vec3 p;
    do {
        p = double(2.0)*Vec3(drand48(), drand48(), drand48()) - Vec3(1, 1, 1);
    } while (p.squared_length() >= 1.0);
    return p;
}*/

class Material {
public:
    virtual bool scatter(const Ray& r_in, const HitRecord& rec, Vec3& attenuation, Ray& scattered) const = 0;
};

class Lambertian: public Material {
public:
    Lambertian(shared_ptr<Texture> a) : albedo(a) {}
    virtual bool scatter(const Ray& r_in, const HitRecord& rec, Vec3& attenuation, Ray& scattered) const {
        Vec3 scatter_direction = rec.normal + random_in_unit_sphere();
        scattered = Ray(rec.p, scatter_direction, double(r_in.get_time()));
        attenuation = albedo->value(rec.u, rec.v, rec.p);
        //attenuation = albedo;
        return true;
    }
public:
    shared_ptr<Texture> albedo;
};

class Metal: public Material {
public:
    Vec3 albedo;
    double fuzz;
    Metal(const Vec3& a, double f) : albedo(a) { if (f < 1) fuzz = f; else fuzz = 1; }
    virtual bool scatter(const Ray& r_in, const HitRecord& rec, Vec3& attenuation, Ray& scattered) const {
        Vec3 reflected = reflect(unit_vector(r_in.get_direction()), rec.normal);
        scattered = Ray(rec.p, reflected + fuzz*random_in_unit_sphere());
        attenuation = albedo;
        return (dot(scattered.get_direction(), rec.normal) > 0);
    }
};

/*class Dielectric: public Material {
public:
    double ref_idx;
    Dielectric(double ri) : ref_idx(ri) {}
    virtual bool scatter(const Ray& r_in, const HitRecord& rec, Vec3& attenuation, Ray& scattered) const {
        Vec3 outward_normal;
        Vec3 reflected = reflect(r_in.get_direction(), rec.normal);
        double ni_over_nt;
        attenuation = Vec3(1.0, 1.0, 1.0);
        Vec3 refracted;
        double reflect_prob;
        double cosine;
        if (dot(r_in.get_direction(), rec.normal) > 0) {
            outward_normal = -rec.normal;
            ni_over_nt = ref_idx;
            cosine = dot(r_in.get_direction(), rec.normal) / r_in.get_direction().length();
            cosine = sqrt(1 - ref_idx*ref_idx*(1 - cosine*cosine));
        }
        else {
            outward_normal = rec.normal;
            ni_over_nt = 1.0 / ref_idx;
            cosine = -dot(r_in.get_direction(), rec.normal) / r_in.get_direction().length();
        }
        if (refract(r_in.get_direction(), outward_normal, ni_over_nt, refracted)) {
            reflect_prob = schlick(cosine, ref_idx);
        }
        else {
            reflect_prob = 1.0;
        }
        if (drand48() < reflect_prob) {
            scattered = Ray(rec.p, reflected);
        }
        else {
            scattered = Ray(rec.p, refracted);
        }
        return true;
    }
};*/

class Dielectric : public Material {
public:
    Dielectric(double ri) : ref_idx(ri) {}

    virtual bool scatter(
            const Ray& r_in, const HitRecord& rec, Vec3& attenuation, Ray& scattered
    ) const {
        attenuation = Vec3(1.0, 1.0, 1.0);
        double etai_over_etat = (rec.front_face) ? (1.0 / ref_idx) : (ref_idx);

        Vec3 unit_direction = unit_vector(r_in.get_direction());
        double cos_theta = ffmin(dot(-unit_direction, rec.normal), 1.0);
        double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
        if (etai_over_etat * sin_theta > 1.0 ) {
            Vec3 reflected = reflect(unit_direction, rec.normal);
            scattered = Ray(rec.p, reflected);
            return true;
        }
        double reflect_prob = schlick(cos_theta, etai_over_etat);
        if (random_double() < reflect_prob)
        {
            Vec3 reflected = reflect(unit_direction, rec.normal);
            scattered = Ray(rec.p, reflected);
            return true;
        }
        Vec3 refracted = refract(unit_direction, rec.normal, etai_over_etat);
        scattered = Ray(rec.p, refracted);
        return true;
    }

public:
    double ref_idx;
};
#endif //RAYTRACER_MATERIAL_H
