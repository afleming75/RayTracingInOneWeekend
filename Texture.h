#ifndef RAYTRACER_TEXTURE_H
#define RAYTRACER_TEXTURE_H
#include "Utilities.h"
#include "Vec3.h"
#include "Perlin.h"

class Texture {
public:
    virtual Vec3 value(double u, double v, const Vec3& p) const = 0;
};

class ConstantTexture : public Texture {
public:
    ConstantTexture() {}
    ConstantTexture(Vec3 c) : color(c) {}

    virtual Vec3 value(double u, double v, const Vec3& p) const {
        return color;
    }
public:
    Vec3 color;
};

class CheckerTexture : public Texture {
public:
    CheckerTexture() {}
    CheckerTexture(shared_ptr<Texture> t0, shared_ptr<Texture> t1): even(t0), odd(t1) {}

    virtual Vec3 value(double u, double v, const Vec3& p) const {
        auto sines = sin(10*p.x())*sin(10*p.y())*sin(10*p.z());
        if (sines < 0)
            return odd->value(u, v, p);
        else
            return even->value(u, v, p);
    }
public:
    shared_ptr<Texture> odd;
    shared_ptr<Texture> even;
};

class NoiseTexture : public Texture {
public:
    NoiseTexture() {}
    NoiseTexture(double sc) : scale(sc) {}
    virtual Vec3 value(double u, double v, const Vec3& p) const {
        return Vec3(1, 1, 1) * 0.5 * (1 + sin(scale*p.z() + 10*noise.turbulence(p)));
    }
    Perlin noise;
    double scale;
};
#endif //RAYTRACER_TEXTURE_H
