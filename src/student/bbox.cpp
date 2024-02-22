
#include "../lib/mathlib.h"
#include "debug.h"

bool BBox::hit(const Ray& ray, Vec2& times) const {

    // TODO (PathTracer):
    // Implement ray - bounding box intersection test
    // If the ray intersected the bounding box within the range given by
    // [times.x,times.y], update times with the new intersection times.

    Vec3 origin = ray.point;
    Vec3 direction = ray.dir;

    float t_min = (min.x - origin.x) / direction.x;
    float t_max = (max.x - origin.x) / direction.x;

    if(t_min > t_max) {
        std::swap(t_min, t_max);
    }

    float ty_min = (min.y - origin.y) / direction.y;
    float ty_max = (max.y - origin.y) / direction.y;

    if(ty_min > ty_max) {
        std::swap(ty_min, ty_max);
    }

    if((t_min > ty_max) || (ty_min > t_max)) {
        return false;
    }

    if(ty_min > t_min) {
        t_min = ty_min;
    }

    if(ty_max < t_max) {
        t_max = ty_max;
    }

    float tz_min = (min.z - origin.z) / direction.z;
    float tz_max = (max.z - origin.z) / direction.z;

    if(tz_min > tz_max) {
        std::swap(tz_min, tz_max);
    }

    if((t_min > tz_max) || (tz_min > t_max)) {
        return false;
    }

    if(tz_min > t_min) {
        t_min = tz_min;
    }

    if(tz_max < t_max) {
        t_max = tz_max;
    }

    times.x = t_min;
    times.y = t_max;
    return true;
}
