
#include "../rays/tri_mesh.h"
#include "debug.h"
#include <iostream>

namespace PT {

BBox Triangle::bbox() const {

    // TODO (PathTracer): Task 2
    // compute the bounding box of the triangle

    // Beware of flat/zero-volume boxes! You may need to
    // account for that here, or later on in BBox::intersect

    BBox box(vertex_list[v0].position, vertex_list[v0].position);
    // box.enclose(vertex_list[v0].position);
    box.enclose(vertex_list[v1].position);
    box.enclose(vertex_list[v2].position);

    return box;
}

Trace Triangle::hit(const Ray& ray) const {

    // Vertices of triangle - has postion and surface normal
    // See rays/tri_mesh.h for a description of this struct
    
    Tri_Mesh_Vert v_0 = vertex_list[v0];
    Tri_Mesh_Vert v_1 = vertex_list[v1];
    Tri_Mesh_Vert v_2 = vertex_list[v2];

    // here just to avoid unused variable warnings, students should remove the following three lines.
    // (void)v_0;
    // (void)v_1;
    // (void)v_2;
    
    // TODO (PathTracer): Task 2
    // Intersect this ray with a triangle defined by the above three points.
    // Intersection should yield a ray t-value, and a hit point (u,v) on the surface of the triangle
    // TODO: Check CCW positions.
    auto det = [] (Vec3 a, Vec3 b, Vec3 c) {
        return dot(cross(a, b), c);
    };

    Vec3 e1 = v_1.position - v_0.position;
    Vec3 e2 = v_2.position - v_0.position;
    Vec3 s = ray.point - v_0.position;
    Vec3 d = ray.dir;
    Trace ret;
    float disc = det(e1, d, e2);
    if (disc == 0) {
        return ret;
    }

    float u = -det(s, e2, d) / disc;
    float v = det(e1, d, s) / disc;
    float t = -det(s, e2, e1) / disc;
    if (u >= 0.0f && u <= 1.0f && v >= 0.0f && v <= 1.0f && u + v <= 1.0f && t >= ray.dist_bounds.x && t <= ray.dist_bounds.y) {
        // You'll need to fill in a "Trace" struct describing information about the hit (or lack of hit)
        ret.origin = ray.point;
        ret.hit = true;     // was there an intersection?
        ret.distance = t;   // at what distance did the intersection occur?
        ret.position = ray.at(t); // where was the intersection?
        ray.dist_bounds.y = t;
        ret.normal = v_0.normal * (1.0f - u - v) + v_1.normal * u + v_2.normal * v;
                            // what was the surface normal at the intersection?
                            // (this should be interpolated between the three vertex normals)
    }
    return ret;
}

Triangle::Triangle(Tri_Mesh_Vert* verts, unsigned int v0, unsigned int v1, unsigned int v2)
    : vertex_list(verts), v0(v0), v1(v1), v2(v2) {
}

void Tri_Mesh::build(const GL::Mesh& mesh) {

    verts.clear();
    triangles.clear();

    for(const auto& v : mesh.verts()) {
        verts.push_back({v.pos, v.norm});
    }

    const auto& idxs = mesh.indices();

    std::vector<Triangle> tris;
    for(size_t i = 0; i < idxs.size(); i += 3) {
        tris.push_back(Triangle(verts.data(), idxs[i], idxs[i + 1], idxs[i + 2]));
    }

    triangles.build(std::move(tris), 4);
}

Tri_Mesh::Tri_Mesh(const GL::Mesh& mesh) {
    build(mesh);
}

Tri_Mesh Tri_Mesh::copy() const {
    Tri_Mesh ret;
    ret.verts = verts;
    ret.triangles = triangles.copy();
    return ret;
}

BBox Tri_Mesh::bbox() const {
    return triangles.bbox();
}

Trace Tri_Mesh::hit(const Ray& ray) const {
    Trace t = triangles.hit(ray);
    return t;
}

size_t Tri_Mesh::visualize(GL::Lines& lines, GL::Lines& active, size_t level,
                           const Mat4& trans) const {
    return triangles.visualize(lines, active, level, trans);
}

} // namespace PT
