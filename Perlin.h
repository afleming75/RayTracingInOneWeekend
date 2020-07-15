#ifndef RAYTRACER_PERLIN_H
#define RAYTRACER_PERLIN_H

#include "Random.h"
#include "Vec3.h"

inline double trilinear_interp(double c[2][2][2], double u, double v, double w) {
    auto accum = 0.0;
    for (int i=0; i < 2; i++)
        for (int j=0; j < 2; j++)
            for (int k=0; k < 2; k++)
                accum += (i*u + (1-i)*(1-u))*
                         (j*v + (1-j)*(1-v))*
                         (k*w + (1-k)*(1-w))*c[i][j][k];

    return accum;
}

inline double perlin_interp(Vec3 c[2][2][2], double u, double v, double w) {
    double uu = u * u*(3 - 2 * u);
    double vv = v * v*(3 - 2 * v);
    double ww = w * w*(3 - 2 * w);
    double accum = 0;
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++) {
                Vec3 weight_v(u - i, v - j, w - k);
                accum += (i*uu + (1 - i)*(1 - uu))*
                         (j*vv + (1 - j)*(1 - vv))*
                         (k*ww + (1 - k)*(1 - ww))*dot(c[i][j][k], weight_v);
            }
    return accum;
}

void permute(int *p, int n) {
    for (int i = n - 1; i > 0; i--) {
        int target = get_random_int(0, i);
        int tmp = p[i];
        p[i] = p[target];
        p[target] = tmp;
    }
    return;
}

static int* perlin_generate_perm() {
    int * p = new int[256];
    for (int i = 0; i < 256; i++)
        p[i] = i;
    permute(p, 256);
    return p;
}

static Vec3* perlin_generate() {
    Vec3 *p = new Vec3[256]; //static const int point_count = 256;
    for (int i = 0; i < 256; ++i) {
        /*double x_random = 2 * (rand() % 100 / float(100)) - 1;
        double y_random = 2 * (rand() % 100 / float(100)) - 1;
        double z_random = 2 * (rand() % 100 / float(100)) - 1;
        p[i] = unit_vector(Vec3(x_random, y_random, z_random));*/
        p[i] = unit_vector(Vec3::random(-1, 1));
    }
    return p;
}

class Perlin {
public:
    /*Perlin() {
        ranvec = perlin_generate();
        perm_x = perlin_generate_perm();
        perm_y = perlin_generate_perm();
        perm_z = perlin_generate_perm();
    }

    ~Perlin() {
        delete[] ranvec;
        delete[] perm_x;
        delete[] perm_y;
        delete[] perm_z;
    }*/

    double noise(const Vec3& p) const {
        double u = p.x() - floor(p.x());
        double v = p.y() - floor(p.y());
        double w = p.z() - floor(p.z());
        int i = floor(p.x());
        int j = floor(p.y());
        int k = floor(p.z());
        Vec3 c[2][2][2];
        for (int di = 0; di < 2; di++)
            for (int dj = 0; dj < 2; dj++)
                for (int dk = 0; dk < 2; dk++)
                    c[di][dj][dk] = ranvec[perm_x[(i + di) & 255] ^ perm_y[(j + dj) & 255] ^ perm_z[(k + dk) & 255]];
        return perlin_interp(c, u, v, w);
    }
    double turbulence(const Vec3& p, int depth = 7) const {
        double accum = 0;
        Vec3 temp_p = p;
        double weight = 1.0;
        for (int i = 0; i < depth; i++) {
            accum += weight * noise(temp_p);
            weight *= 0.5;
            temp_p *= 2;
        }
        return fabs(accum);
    }

    static Vec3 *ranvec;
    static int *perm_x;
    static int *perm_y;
    static int *perm_z;
};

Vec3 *Perlin::ranvec = perlin_generate();
int *Perlin::perm_x = perlin_generate_perm();
int *Perlin::perm_y = perlin_generate_perm();
int *Perlin::perm_z = perlin_generate_perm();

#endif //RAYTRACER_PERLIN_H
