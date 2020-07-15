cmake_minimum_required(VERSION 3.16)
project(RayTracer)

set(CMAKE_CXX_STANDARD 14)

add_executable(RayTracer main.cpp HittableList.h Material.h MovingSphere.h Random.h Ray.h Sphere.h Vec3.h Aabb.h BvhNode.h Camera.h Hittable.h Utilities.h Texture.h Perlin.h)