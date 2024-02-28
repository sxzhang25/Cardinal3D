
#include "../rays/shapes.h"
#include "debug.h"

namespace PT {

const char* Shape_Type_Names[(int)Shape_Type::count] = {"None", "Sphere"};

BBox Sphere::bbox() const {

    BBox box;
    box.enclose(Vec3(-radius));
    box.enclose(Vec3(radius));
    return box;
}

Trace Sphere::hit(const Ray& ray) const {

    // TODO (PathTracer): Task 2
    // Intersect this ray with a sphere of radius Sphere::radius centered at the origin.

    // If the ray intersects the sphere twice, ret should
    // represent the first intersection, but remember to respect
    // ray.dist_bounds! For example, if there are two intersections,
    // but only the _later_ one is within ray.dist_bounds, you should
    // return that one!

    Trace ret;
    ret.origin = ray.point;

    Vec3 o = ray.point;
    Vec3 d = ray.dir;
    float b = dot(o, d);
    float det = b * b - dot(o, o) + radius * radius;
    if (det < 0.0f) return ret;

    float t1 = -b - sqrt(det);
    float t2 = -b + sqrt(det);
    if (t1 >= ray.dist_bounds.x && t1 <= ray.dist_bounds.y) {
        ret.hit = true;
        ret.position = ray.at(t1);
        ret.distance = t1;
        ret.normal = ret.position.unit();
        ray.dist_bounds.y = t1;
        return ret;
    }
    if (t2 >= ray.dist_bounds.x && t2 <= ray.dist_bounds.y) {
        ret.hit = true;
        ret.position = ray.at(t2);
        ret.distance = t2;
        ret.normal = ret.position.unit();
        ray.dist_bounds.y = t2;
        return ret;
    }
    return ret;
}

} // namespace PT
