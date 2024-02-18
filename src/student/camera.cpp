
#include "../util/camera.h"
#include "../rays/samplers.h"
#include "debug.h"

Ray Camera::generate_ray(Vec2 screen_coord) const {

    // TODO (PathTracer): Task 1
    //
    // The input screen_coord is a normalized screen coordinate [0,1]^2
    //
    // You need to transform this 2D point into a 3D position on the sensor plane, which is
    // located one unit away from the pinhole in camera space (aka view space).
    //
    // You'll need to compute this position based on the vertial field of view
    // (vert_fov) of the camera, and the aspect ratio of the output image (aspect_ratio).
    //
    // Tip: compute the ray direction in view space and use
    // the camera space to world space transform (iview) to transform the ray back into world space.

    float vs_height = 2.0f * tan(0.5f * M_PI * vert_fov / 180.0f);
    float vs_width = aspect_ratio * vs_height;
    Vec3 vs_dir = Vec3(vs_width * (screen_coord.x - 0.5f), vs_height * (screen_coord.y - 0.5f), -1.0f);
    Vec3 cs_dir = iview * vs_dir;
    Ray r = Ray(position, cs_dir);

    return r;
}
