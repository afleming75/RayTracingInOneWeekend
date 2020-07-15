
#ifndef RAYTRACER_CAMERA_H
#define RAYTRACER_CAMERA_H
#include "Random.h"
#include "Ray.h"

Vec3 random_in_unit_disk() {
    while (true) {
        auto p = Vec3(get_random_double(-1, 1), get_random_double(-1, 1), 0);
        if(p.squared_length() >= 1) continue;
        return p;
    }
    /*Vec3 p;
    do {
        p = 2.0*Vec3(drand48(), drand48(), 0) - Vec3(1, 1, 0);
    } while (dot(p, p) >= 1.0);
    return p;*/
}

class Camera {
public:
    Vec3 lower_left_corner;
    Vec3 horizontal;
    Vec3 vertical;
    Vec3 origin;
    Vec3 u, v, w;
    double time0, time1; //shutter open and close times
    double lens_radius;
    Camera(Vec3 lookfrom, Vec3 lookat, Vec3 vup, double vfov, double aspect, double aperture, double focus_dist,
           double t0 = 0, double t1 = 0) {
        //vfov - top to bottom, in degrees
        lens_radius = aperture / 2;
        time0 = t0;
        time1 = t1;
        auto theta = vfov*M_PI/180; //degrees to radians
        auto half_height = tan(theta/2);
        auto half_width = aspect*half_height;
        origin = lookfrom;
        w = unit_vector(lookfrom - lookat);
        u = unit_vector(cross(vup, w));
        v = cross(w, u);
        lower_left_corner = origin - half_width*focus_dist*u - half_height*focus_dist*v - focus_dist*w;
        horizontal = 2*half_width*focus_dist*u;
        vertical = 2*half_height*focus_dist*v;
    }
    Ray get_ray(double s, double t) {
        Vec3 rd = lens_radius*random_in_unit_disk();
        Vec3 offset = u*rd.x() + v*rd.y();
        return Ray(origin + offset, lower_left_corner + s*horizontal + t*vertical - origin - offset,
                   get_random_double(time0, time1));
    }

};
#endif //RAYTRACER_CAMERA_H
