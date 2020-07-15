#ifndef RAYTRACER_VEC3_H
#define RAYTRACER_VEC3_H
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "Random.h"
#include "Utilities.h"

class Vec3 {
public:
    double e[3];
    Vec3() {}
    Vec3(double e0, double e1, double e2) { e[0] = e0;	e[1] = e1; e[2] = e2; }
    inline double x() const { return e[0]; }
    inline double y() const { return e[1]; }
    inline double z() const { return e[2]; }
    inline double r() const { return e[0]; }
    inline double g() const { return e[1]; }
    inline double b() const { return e[2]; }

    inline const Vec3& operator+() const { return *this; } //get copy as reference
    inline Vec3 operator-() const { return Vec3(-e[0], -e[1], -e[2]); }
    inline double operator[](int i) const { return e[i]; }
    inline double& operator[](int i) { return e[i]; };

    inline Vec3& operator+=(const Vec3 &v2);
    inline Vec3& operator-=(const Vec3 &v2);
    inline Vec3& operator*=(const Vec3 &v2);
    inline Vec3& operator/=(const Vec3 &v2);
    inline Vec3& operator*=(const double t);
    inline Vec3& operator/=(const double t);

    inline double length() const {
        return sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]); }
    inline double squared_length() const {
        return e[0]*e[0] + e[1]*e[1] + e[2]*e[2]; }
    inline void make_unit_vector();
    void write_color(std::ofstream &out, int samples_per_pixel);

    inline static Vec3 random() {
        return Vec3(random_double(), random_double(), random_double());
    }

    inline static Vec3 random(double min, double max) {
        return Vec3(get_random_double(min,max), get_random_double(min,max), get_random_double(min,max));
    }
};

inline std::istream& operator>>(std::istream &is, Vec3 &t) {
    is >> t.e[0] >> t.e[1] >> t.e[2];
    return is;
}

inline std::ostream& operator<<(std::ostream &os, const Vec3 &t) {
    os << t.e[0] << " " << t.e[1] << " " << t.e[2];
    return os;
}

inline void Vec3::make_unit_vector() {
    double k = 1.0 / sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
    e[0] *= k;
    e[1] *= k;
    e[2] *= k;
}

inline Vec3 operator+(const Vec3 &v1, const Vec3 &v2) {
    return Vec3(v1.e[0] + v2.e[0], v1.e[1] + v2.e[1], v1.e[2] + v2.e[2]);
}

inline Vec3 operator-(const Vec3 &v1, const Vec3 &v2) {
    return Vec3(v1.e[0] - v2.e[0], v1.e[1] - v2.e[1], v1.e[2] - v2.e[2]);
}

inline Vec3 operator*(const Vec3 &v1, const Vec3 &v2) {
    return Vec3(v1.e[0] * v2.e[0], v1.e[1] * v2.e[1], v1.e[2] * v2.e[2]);
}

inline Vec3 operator/(const Vec3 &v1, const Vec3 &v2) {
    return Vec3(v1.e[0] / v2.e[0], v1.e[1] / v2.e[1], v1.e[2] / v2.e[2]);
}

inline Vec3 operator*(double t, const Vec3 &v) {
    return Vec3(t*v.e[0], t*v.e[1], t*v.e[2]);
}


inline Vec3 operator/(Vec3 v, double t) {
    return Vec3(v.e[0]/t, v.e[1]/t, v.e[2]/t);
}

inline Vec3 operator*(const Vec3 &v, double t) {
    return Vec3(t*v.e[0], t*v.e[1], t*v.e[2]);
}

inline double dot(const Vec3 &v1, const Vec3 &v2) {
    return v1.e[0]*v2.e[0] + v1.e[1]*v2.e[1] + v1.e[2]*v2.e[2];
}

inline Vec3 cross(const Vec3 &v1, const Vec3 &v2) {
    return Vec3( v1.e[1]*v2.e[2] - v1.e[2]*v2.e[1],
                 v1.e[2]*v2.e[0] - v1.e[0]*v2.e[2],
                 v1.e[0]*v2.e[1] - v1.e[1]*v2.e[0]);
}

inline Vec3& Vec3::operator+=(const Vec3 &v) {
    e[0] += v.e[0];
    e[1] += v.e[1];
    e[2] += v.e[2];
    return *this;
}

inline Vec3& Vec3::operator*=(const Vec3 &v) {
    e[0] *= v.e[0];
    e[1] *= v.e[1];
    e[2] *= v.e[2];
    return *this;
}

inline Vec3& Vec3::operator/=(const Vec3 &v) {
    e[0] /= v.e[0];
    e[1] /= v.e[1];
    e[2] /= v.e[2];
    return *this;
}

inline Vec3& Vec3::operator-=(const Vec3 &v) {
    e[0] -= v.e[0];
    e[1] -= v.e[1];
    e[2] -= v.e[2];
    return *this;
}

inline Vec3& Vec3::operator*=(const double t) {
    e[0] *= t;
    e[1] *= t;
    e[2] *= t;
    return *this;
}

inline Vec3& Vec3::operator/=(const double t) {
    double k = 1.0/t;
    e[0] *= k;
    e[1] *= k;
    e[2] *= k;
    return *this;
}

inline Vec3 random_unit_vector() {
    Vec3 p;
    p = Vec3(random_double(), random_double(), random_double());
    return p / p.length();
}

Vec3 random_in_unit_sphere() {
    while (true) {
        auto p = Vec3::random(-1,1);
        if (p.squared_length() >= 1) continue;
        return p;
    }
}

Vec3 random_in_hemisphere(const Vec3& normal) {
    Vec3 in_unit_sphere = random_in_unit_sphere();
    if (dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}

inline Vec3 unit_vector(Vec3 v) {
    return v / v.length();
}

void Vec3::write_color(std::ofstream &out, int samples_per_pixel) {
    // Divide the color total by the number of samples and gamma-correct
    // for a gamma value of 2.0.
    auto scale = 1.0 / samples_per_pixel;
    auto r = sqrt(scale * e[0]);
    auto g = sqrt(scale * e[1]);
    auto b = sqrt(scale * e[2]);

    // Write the translated [0,255] value of each color component.
    out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
}
#endif //RAYTRACER_VEC3_H
