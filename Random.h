
#ifndef RAYTRACER_RANDOM_H
#define RAYTRACER_RANDOM_H
#include <cstdlib>
#include <random>
#include <functional>
/*inline double random_double() {
    return rand() / (RAND_MAX + 1.0);
}

inline double get_random_double(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}*/
inline double get_random_double(double min, double max) {
    static std::uniform_real_distribution<double> distribution(min, max);
    static std::mt19937 generator;
    static std::function<double()> rand_generator =
            std::bind(distribution, generator);
    return rand_generator();
}

inline double random_double() {
    return get_random_double(0.0, 1.0);
}

inline double get_random_int(int fMin, int fMax) {
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_int_distribution<int>  distr(fMin, fMax);

    return distr(generator);
}

#endif //RAYTRACER_RANDOM_H
