
#ifndef RAYTRACER_HITTABLELIST_H
#define RAYTRACER_HITTABLELIST_H
#include "Hittable.h"
#include <algorithm>
#include <vector>
#include <memory>

class HittableList: public Hittable {
public:

    HittableList() {}
    HittableList(std::shared_ptr<Hittable> object) { add(object); }

    void clear() { objects.clear(); }
    void add(std::shared_ptr<Hittable> object) { objects.push_back(object); }

    virtual bool hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const;
    virtual bool bounding_box(double t0, double t1, Aabb& output_box) const;

public:
    std::vector<std::shared_ptr<Hittable>> objects;

};


bool HittableList::bounding_box(double t0, double t1, Aabb& output_box) const {
    if (objects.empty())
        return false;
    Aabb temp_box;
    bool first_box = true;
    for (const auto& object : objects) {
        if (!object->bounding_box(t0, t1, temp_box))
            return false;
        output_box = first_box ? temp_box : surrounding_box(output_box, temp_box);
        first_box = false;
    }
    return true;
}

bool HittableList::hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const {
    HitRecord temp_rec;
    bool hit_anything = false;
    double closest_so_far = t_max;
    for (const auto& object : objects) {
            if (object->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }
    return hit_anything;
}
#endif //RAYTRACER_HITTABLELIST_H

/*#ifndef RAYTRACER_HITTABLELIST_H
#define RAYTRACER_HITTABLELIST_H
#include "Hittable.h"
#include <algorithm>

class HittableList: public Hittable {
public:
    HittableList();
    int list_size;
    Hittable **list;
    HittableList(Hittable **l, int n) {list = l; list_size = n; }
    virtual bool hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const;
    virtual bool bounding_box(double t0, double t1, Aabb& output_box) const;
    //bool empty() const;
    //HittableList add();
};


bool HittableList::bounding_box(double t0, double t1, Aabb& output_box) const {
    if (list_size < 1) return false; //(objects.empty())

    Aabb temp_box;
    bool first_box = list[0]->bounding_box(t0, t1, temp_box);
    if(!first_box)
        return false;
    else
        output_box = temp_box;
    for (int i = 1; i < list_size; i++) {
        if (list[i]->bounding_box(t0, t1, temp_box)) {
            output_box = surrounding_box(output_box, temp_box);
        } else
            return false;
    }
    return true;
}

bool HittableList::hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const {
    HitRecord temp_rec;
    bool hit_anything = false;
    double closest_so_far = t_max;
    for (int i = 0; i < list_size; i++) {
        if (list[i]->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }
    return hit_anything;
}
#endif //RAYTRACER_HITTABLELIST_H */
