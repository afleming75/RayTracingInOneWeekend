#include <iostream>
#include <fstream>
#include "float.h"
#include <memory>
#include "Camera.h"
#include "HittableList.h"
#include "Material.h"
#include "Sphere.h"
#include "MovingSphere.h"
#include "Vec3.h"

using std::shared_ptr;
using std::make_shared;



/*Hittable *random_scene() {
    int n = 500;
    Hittable **list = make_shared<Hittable>[n + 1];
    list[0] = new Sphere(Vec3(0, -1000, 0), 1000, new Lambertian(Vec3(0.5, 0.5, 0.5)));
    int i = 1;
    for (int a = -10; a < 10; a++) {
        for(int b = -10; b < 10; b++) {
            auto choose_mat = random_double();
            Vec3 center(a + 0.9*random_double(), 0.2, b + 0.9*random_double());
            if ((center - Vec3(4, 0.2, 0)).length() > 0.9) {
                if(choose_mat < 0.8) { //diffuse
                    //Vec3 temp = new Vec3(0.0, 0.0, 0.0);
                    auto albedo = (random_unit_vector()*random_unit_vector());
                    list[i++] = make_shared<MovingSphere>(center,
                                                 center + Vec3(0, get_random_double(0,  0.5), 0), 0.0, 1.0, 0.2,
                                                 new Lambertian(albedo));
                } else if (choose_mat < 0.95) { //metal
                    auto albedo = random_unit_vector();
                    auto fuzz = get_random_double(0.0, 0.5);
                    list[i++] = new Sphere(center, 0.2, new Metal(albedo, fuzz));
                } else { //glass
                    list[i++] = new Sphere(center, 0.2, new Dielectric(1.5));
                }
            }
        }
    }
    list[i++] = new Sphere(Vec3(0, 1, 0), 1.0, new Dielectric(1.5));
    list[i++] = new Sphere(Vec3(-4, 1, 0), 1.0, new Lambertian(Vec3(0.4, 0.2, 0.1)));
    list[i++] = new Sphere(Vec3(4, 1, 0), 1.0, new Metal(Vec3(0.7, 0.6, 0.5), 0.0));

    return new HittableList(list, i);
}*/

/*Hittable *random_scene() {
	int n = 500;
	Hittable **list = new Hittable*[n + 1];
	list[0] = new Sphere(Vec3(0, -1000, 0), 1000, new Lambertian(Vec3(0.5, 0.5, 0.5)));
	int i = 1;
	for(int a = -11; a < 11; a++) {
		for(int b = -11; b < 11; b++) {
			double choose_mat = drand48();
			Vec3 center(a + 0.9*drand48(), 0.2, b + 0.9*drand48());
			if ((center - Vec3(4, 0.2, 0)).length() > 0.9) {
				if (choose_mat < 0.8) { //diffuse
					list[i++] = new Sphere(center, 0.2, new Lambertian(Vec3(drand48()*drand48(), drand48()*drand48(), drand48()*drand48())));
				}
				else if(choose_mat < 0.95) { //metal
					list[i++] = new Sphere(center, 0.2, new Metal(Vec3(0.5*(1 + drand48()), 0.5*(1 + drand48()), 0.5*(1 + drand48())), 0.5*drand48()));
				}
				else { //glass
					list[i++] = new Sphere(center, 0.2, new Dielectric(1.5));
				}
			}
		}
	}
	list[i++] = new Sphere(Vec3(0, 1, 0), 1.0, new Dielectric(1.5));
	list[i++] = new Sphere(Vec3(-4, 1, 0), 1.0, new Lambertian(Vec3(0.4, 0.2, 0.1)));
	list[i++] = new Sphere(Vec3(4, 1, 0), 1.0, new Metal(Vec3(0.7, 0.6, 0.5), 0.0));

	return world;
}*/

HittableList random_scene() {
    HittableList world;

    world.add(make_shared<Sphere>(
            Vec3(0,-1000,0), 1000, make_shared<Lambertian>(make_shared<ConstantTexture>(Vec3(0.5, 0.5, 0.5)))));

    int i = 1;
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_double();
            Vec3 center(a + 0.9*random_double(), 0.2, b + 0.9*random_double());
            if ((center - Vec3(4, 0.2, 0)).length() > 0.9) {
                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = random_unit_vector() * random_unit_vector();
                    world.add(
                            make_shared<Sphere>(center, 0.2, make_shared<Lambertian>(make_shared<ConstantTexture>(albedo))));
                } else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = random_unit_vector();//(.5, 1);
                    auto fuzz = get_random_double(0, .5);
                    world.add(
                            make_shared<Sphere>(center, 0.2, make_shared<Metal>(albedo, fuzz)));
                } else {
                    // glass
                    world.add(make_shared<Sphere>(center, 0.2, make_shared<Dielectric>(1.5)));
                }
            }
        }
    }

    world.add(make_shared<Sphere>(Vec3(0, 1, 0), 1.0, make_shared<Dielectric>(1.5)));

    world.add(
            make_shared<Sphere>(Vec3(-4, 1, 0), 1.0, make_shared<Lambertian>(make_shared<ConstantTexture>(Vec3(0.4, 0.2, 0.1)))));

    world.add(
            make_shared<Sphere>(Vec3(4, 1, 0), 1.0, make_shared<Metal>(Vec3(0.7, 0.6, 0.5), 0.0)));

    return world;
}

HittableList three_spheres() {
    HittableList world;
    world.add(make_shared<Sphere>(Vec3(0, 1, 0), 1.0, make_shared<Dielectric>(1.5)));
    world.add(make_shared<Sphere>(Vec3(-4, 1, 0), 1.0,
                                  make_shared<Lambertian>(make_shared<ConstantTexture>(Vec3(0.4, 0.2, 0.1)))));
    world.add(make_shared<Sphere>(Vec3(4, 1, 0), 1.0, make_shared<Metal>(Vec3(0.7, 0.6, 0.5), 0.0)));
    return world;
}

HittableList random_moving_scene() {
    HittableList world;

    auto checker = make_shared<CheckerTexture>(
            make_shared<ConstantTexture>(Vec3(0.2, 0.3, 0.1)),
            make_shared<ConstantTexture>(Vec3(0.9, 0.9, 0.9))
    );
    world.add(make_shared<Sphere>(
            Vec3(0,-1000,0), 1000, make_shared<Lambertian>(checker)));

    int i = 1;
    for (int a = -10; a < 10; a++) {
        for (int b = -10; b < 10; b++) {
            auto choose_mat = random_double();
            Vec3 center(a + 0.9*random_double(), 0.2, b + 0.9*random_double());
            if ((center - Vec3(4, .2, 0)).length() > 0.9) {
                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = Vec3::random() * Vec3::random();
                    world.add(make_shared<MovingSphere>(
                            center, center + Vec3(0, get_random_double(0,.5), 0), 0.0, 1.0, 0.2,
                            make_shared<Lambertian>(make_shared<ConstantTexture>(albedo))));
                } else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = Vec3::random(.5, 1);
                    auto fuzz = get_random_double(0, .5);
                    world.add(
                            make_shared<Sphere>(center, 0.2, make_shared<Metal>(albedo, fuzz)));
                } else {
                    // glass
                    world.add(make_shared<Sphere>(center, 0.2, make_shared<Dielectric>(1.5)));
                }
            }
        }
    }

    world.add(make_shared<Sphere>(Vec3(0, 1, 0), 1.0, make_shared<Dielectric>(1.5)));
    world.add(make_shared<Sphere>(
            Vec3(-4, 1, 0), 1.0, make_shared<Lambertian>(make_shared<ConstantTexture>(Vec3(0.4, 0.2, 0.1)))));
    world.add(make_shared<Sphere>(
            Vec3(4, 1, 0), 1.0, make_shared<Metal>(Vec3(0.7, 0.6, 0.5), 0.0)));

    return world;
}

HittableList two_perlin_spheres() {
    HittableList objects;

    auto pertext = make_shared<NoiseTexture>(4);//4
    objects.add(make_shared<Sphere>(Vec3(0,-1000, 0), 1000, make_shared<Lambertian>(pertext)));
    objects.add(make_shared<Sphere>(Vec3(0, 2, 0), 2, make_shared<Lambertian>(pertext)));

    return objects;
}
/*Vec3 color(const Ray& r, const Hittable& world, int depth) {
    HitRecord rec;
    if (world.hit(r, 0.001, MAXFLOAT, rec)) {
        Ray scattered;
        Vec3 attenuation;
        if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
            return attenuation*color(scattered, world, depth + 1);
        }
        else {
            return Vec3(0, 0, 0);
        }
    }
    else {
        Vec3 unit_direction = unit_vector(r.get_direction());
        double t = 0.5*(unit_direction.y() + 1.0);
        return (1.0 - t)*Vec3(1.0, 1.0, 1.0) + t*Vec3(0.5, 0.7, 1.0);
    }
}*/

Vec3 color(const Ray& r, const Hittable& world, int depth) {
    HitRecord rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return Vec3(0,0,0);

    if (world.hit(r, 0.001, infinity(), rec)) {
        Ray scattered;
        Vec3 attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * color(scattered, world, depth-1);
        return Vec3(0,0,0);
    }

    Vec3 unit_direction = unit_vector(r.get_direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*Vec3(1.0, 1.0, 1.0) + t*Vec3(0.5, 0.7, 1.0);
}

/*Vec3 color(const Ray& r, const Hittable& world, int depth) {
    HitRecord rec;
    if (depth <= 0)
        return Vec3(0,0,0);
    if (world.hit(r, 0.001, infinity(), rec)) {
        Vec3 target = rec.p + rec.normal + random_in_hemisphere(rec.normal);
        return 0.5 * color(Ray(rec.p, target - rec.p), world, depth-1);
    }

    Vec3 unit_direction = unit_vector(r.get_direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0 - t)*Vec3(1.0, 1.0, 1.0) + t*Vec3(0.5, 0.7, 1.0);
}*/

int main() {
    std::ofstream imageFile;
    imageFile.open ("C:\\Users\\ChangeNameLater\\Desktop\\finalImage.ppm");
    int image_width = 1200;
    int image_height = 800;
    int samplesPerPixel = 10;
    const int max_depth = 50;
    imageFile << "P3\n" << image_width << " " << image_height << "\n255\n";
    HittableList world = two_perlin_spheres();
    //double R = cos(M_PI / 4);
    const auto aspect_ratio = double(image_width) / image_height;
    Vec3 lookfrom(13, 2, 3);
    Vec3 lookat(0, 0, 0);
    Vec3 vup(0, 1, 0);
    double dist_to_focus = 10.0; //(lookfrom - lookat).length();
    double aperture = 0.0;

    Camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);

    for (int j = image_height-1; j >= 0; j--){
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; i++){
            Vec3 col(0, 0, 0);
            for (int s = 0; s < samplesPerPixel; s++) {
                double u = double(i + random_double()) / double(image_width);
                double v = double(j + random_double()) / double(image_height);
                Ray r = cam.get_ray(u, v);
                col += color(r, world, max_depth);
            }
            col.write_color(imageFile, samplesPerPixel);
            /*col /= double(samplesPerPixel);
            col = Vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
            int ir = int(255.99*col[0]);
            int ig = int(255.99*col[1]);
            int ib = int(255.99*col[2]);
            imageFile << ir << " " << ig << " " << ib << "\n";*/
        }
    }
    std::cerr << "\nDone.\n";
    imageFile.close();
}

