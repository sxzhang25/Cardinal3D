
#include "../rays/bvh.h"
#include "debug.h"
#include <stack>

namespace PT {

// helper function to compute the bucket index
template<typename Primitive>
size_t BVH<Primitive>::compute_bucket(BBox box, Vec3 center, int axis) {
    float pos = (center[axis] - box.min[axis]) / (box.max[axis] - box.min[axis]);
    size_t idx = (size_t)(pos * (float)PARTS);
    return std::min<size_t>(idx, PARTS - 1);
}

// construct BVH hierarchy given a vector of prims
template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size) {

    // NOTE (PathTracer):
    // This BVH is parameterized on the type of the primitive it contains. This allows
    // us to build a BVH over any type that defines a certain interface. Specifically,
    // we use this to both build a BVH over triangles within each Tri_Mesh, and over
    // a variety of Objects (which might be Tri_Meshes, Spheres, etc.) in Pathtracer.
    //
    // The Primitive interface must implement these two functions:
    //      BBox bbox() const;
    //      Trace hit(const Ray& ray) const;
    // Hence, you may call bbox() and hit() on any value of type Primitive.

    // Keep these two lines of code in your solution. They clear the list of nodes and
    // initialize member variable 'primitives' as a vector of the scene prims
    nodes.clear();
    primitives = std::move(prims);

    // TODO (PathTracer): Task 3
    // Modify the code ahead to construct a BVH from the given vector of primitives and maximum leaf
    // size configuration.
    //
    // Please use the SAH as described in class.  We recomment the binned build from lecture.
    // In general, here is a rough sketch:
    //
    //  For each axis X,Y,Z:
    //     Try possible splits along axis, evaluate SAH for each
    //  Take minimum cost across all axes.
    //  Partition primitives into a left and right child group
    //  Compute left and right child bboxes
    //  Make the left and right child nodes.
    //
    //
    // While a BVH is conceptually a tree structure, the BVH class uses a single vector (nodes)
    // to store all the nodes. Therefore, BVH nodes don't contain pointers to child nodes,
    // but rather the indices of the
    // child nodes in this array. Hence, to get the child of a node, you have to
    // look up the child index in this vector (e.g. nodes[node.l]). Similarly,
    // to create a new node, don't allocate one yourself - use BVH::new_node, which
    // returns the index of a newly added node.
    //
    // As an example of how to make nodes, the starter code below builds a BVH with a
    // root node that encloses all the primitives and its two descendants at Level 2.
    // For now, the split is hardcoded such that the first primitive is put in the left
    // child of the root, and all the other primitives are in the right child.
    // There are no further descendants.

    // edge case
    if(primitives.empty()) {
        return;
    }

    typedef struct {
        BBox bbox;
        size_t prim_count;
        float total_area;
    } Bucket;

    // Replace these
    std::stack<size_t> S;
    BBox box;
    for(const Primitive& prim : primitives) box.enclose(prim.bbox());

    root_idx = new_node(box, 0, primitives.size(), 0, 0);
    S.push(root_idx);

    while(!S.empty()) {
        size_t curIdx = S.top();
        S.pop();
        Node& curNode = nodes[curIdx];
        size_t pstart = curNode.start;
        size_t pend = pstart + curNode.size;
        if(curNode.size <= max_leaf_size) {
            continue;
        }
        std::vector<Bucket> buckets(PARTS);
        BBox minBoxA;
        BBox minBoxB;
        int minAxis = 0;
        size_t minACount = 0;
        size_t minBCount = 0;
        int minCut = 0;
        BBox curBox = curNode.bbox;
        float minSAH = FLT_MAX;
        for(int axis = 0; axis < 3; axis++) {
            for(int i = 0; i < PARTS; i++) {
                buckets[i].bbox.reset();
                buckets[i].prim_count = 0;
                buckets[i].total_area = 0.0f;
            }
            for(size_t i = pstart; i < pend; i++) {
                Primitive& prim = primitives.at(i);
                size_t bid = compute_bucket(curBox, prim.bbox().center(), axis);
                Bucket& B = buckets.at(bid);
                B.bbox.enclose(prim.bbox());
                B.prim_count++;
            }
            for(size_t i = 1; i < PARTS; i++) {
                BBox boxA = BBox();
                BBox boxB = BBox();
                size_t ACount = 0;
                size_t BCount = 0;
                for(size_t j = 0; j < i; j++) {
                    boxA.enclose(buckets.at(j).bbox);
                    ACount += buckets.at(j).prim_count;
                }
                for(size_t j = i; j < PARTS; j++) {
                    boxB.enclose(buckets.at(j).bbox);
                    BCount += buckets.at(j).prim_count;
                }
                float SAH =
                    (float)ACount * boxA.surface_area() + (float)BCount * boxB.surface_area();
                if(SAH < minSAH) {
                    // save some info
                    minSAH = SAH;
                    minBoxA = boxA;
                    minBoxB = boxB;
                    minBCount = BCount;
                    minACount = ACount;
                    minCut = (int)i;
                    minAxis = axis;
                }
            }
        }
        BBox boxA = minBoxA;
        BBox boxB = minBoxB;
        size_t ACount = minACount;
        size_t BCount = minBCount;
        auto part = [curBox, minAxis, minCut](Primitive& prim) {
            float pos = (prim.bbox().center()[minAxis] - curBox.min[minAxis]) /
                        (curBox.max[minAxis] - curBox.min[minAxis]);
            int idx = (size_t)(pos * (float)PARTS);
            return idx < minCut;
        };
        std::partition(primitives.begin() + pstart, primitives.begin() + pend, part);
        size_t lchild = new_node(boxA, pstart, ACount, 0, 0);
        size_t rchild = new_node(boxB, pstart + ACount, BCount, 0, 0);
        nodes[curIdx].l = lchild;
        nodes[curIdx].r = rchild;
        S.push(lchild);
        S.push(rchild);
    }
}

template<typename Primitive>
Trace BVH<Primitive>::hit(const Ray& ray) const {

    // TODO (PathTracer): Task 3
    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.

    // implementing the front-to-back traversal described in lecture
    static size_t curr_node = 0;
    Trace closest;

    if(nodes[curr_node].is_leaf()) {
        size_t start = nodes[curr_node].start;
        size_t size = nodes[curr_node].size;
        for(size_t i = start; i < start + size; i++) {
            closest = Trace::min(closest, primitives[i].hit(ray));
        }

    } else {
        
        // left and right child indices
        size_t left_node = nodes[curr_node].l;
        size_t right_node = nodes[curr_node].r;

        // left and right child hit
        Vec2 t1, t2, t_second;
        primitives[left_node].bbox().hit(ray, t1);
        primitives[right_node].bbox().hit(ray, t2);

        // indices of the closest and second closest children
        size_t first, second;
        first = left_node;
        second = right_node;
        t_second = t2;

        // swap if t1 is greater
        if(t1.x >= t2.x) {
            first = right_node;
            second = left_node;
            t_second = t1;
        }

        curr_node = first;
        closest = hit(ray);

        // if the second closest is closer than the closest, traverse to the second closest
        Trace second_closest_hit;
        if(t_second.x < closest.distance) {
            curr_node = second;
            second_closest_hit = hit(ray);
        }

        if(second_closest_hit.hit) {
            closest = second_closest_hit;
        }
    }

    return closest;
}


template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
    build(std::move(prims), max_leaf_size);
}

template<typename Primitive>
BVH<Primitive> BVH<Primitive>::copy() const {
    BVH<Primitive> ret;
    ret.nodes = nodes;
    ret.primitives = primitives;
    ret.root_idx = root_idx;
    return ret;
}

template<typename Primitive>
bool BVH<Primitive>::Node::is_leaf() const {
    return l == r;
}

template<typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
    Node n;
    n.bbox = box;
    n.start = start;
    n.size = size;
    n.l = l;
    n.r = r;
    nodes.push_back(n);
    return nodes.size() - 1;
}

template<typename Primitive>
BBox BVH<Primitive>::bbox() const {
    return nodes[root_idx].bbox;
}

template<typename Primitive>
std::vector<Primitive> BVH<Primitive>::destructure() {
    nodes.clear();
    return std::move(primitives);
}

template<typename Primitive>
void BVH<Primitive>::clear() {
    nodes.clear();
    primitives.clear();
}

template<typename Primitive>
size_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, size_t level,
                                 const Mat4& trans) const {

    std::stack<std::pair<size_t, size_t>> tstack;
    tstack.push({root_idx, 0});
    size_t max_level = 0;

    if(nodes.empty()) return max_level;

    while(!tstack.empty()) {

        auto [idx, lvl] = tstack.top();
        max_level = std::max(max_level, lvl);
        const Node& node = nodes[idx];
        tstack.pop();

        Vec3 color = lvl == level ? Vec3(1.0f, 0.0f, 0.0f) : Vec3(1.0f);
        GL::Lines& add = lvl == level ? active : lines;

        BBox box = node.bbox;
        box.transform(trans);
        Vec3 min = box.min, max = box.max;

        auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

        edge(min, Vec3{max.x, min.y, min.z});
        edge(min, Vec3{min.x, max.y, min.z});
        edge(min, Vec3{min.x, min.y, max.z});
        edge(max, Vec3{min.x, max.y, max.z});
        edge(max, Vec3{max.x, min.y, max.z});
        edge(max, Vec3{max.x, max.y, min.z});
        edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
        edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
        edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
        edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
        edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
        edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

        if(node.l && node.r) {
            tstack.push({node.l, lvl + 1});
            tstack.push({node.r, lvl + 1});
        } else {
            for(size_t i = node.start; i < node.start + node.size; i++) {
                size_t c = primitives[i].visualize(lines, active, level - lvl, trans);
                max_level = std::max(c, max_level);
            }
        }
    }
    return max_level;
}

} // namespace PT
