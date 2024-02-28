
#include "../rays/bvh.h"
#include "debug.h"
#include <stack>

namespace PT {

// helper function to compute the bucket index
template<typename Primitive>
size_t BVH<Primitive>::compute_bucket(BBox box, Vec3 center, int axis, int num_buckets) {
    float length = box.max[axis] - box.min[axis];
    float pos = center[axis] - box.min[axis];
    float proportion = pos / length;
    size_t idx = (size_t)(proportion * (float)num_buckets);
    size_t max_idx = (size_t)num_buckets - 1;

    if(idx > max_idx) {
        return max_idx;
    }
    return idx;
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
        size_t num_prim;
        float total_area;
    } Bucket;

    // number of buckets is usally no more than 32
    const int num_buckets = 32;

    std::stack<size_t> bbox_stack;
    BBox root_bbox;

    // initialize the bounding box to enclose all the primitives
    for (size_t i = 0; i < primitives.size(); i++) {
        root_bbox.enclose(primitives[i].bbox());
    }

    root_idx = new_node(root_bbox, 0, primitives.size(), 0, 0);
    bbox_stack.push(root_idx);

    while(!bbox_stack.empty()) {
        size_t curr_id = bbox_stack.top();
        bbox_stack.pop();
        Node& curr_node = nodes[curr_id];

        size_t prim_start = curr_node.start;
        size_t prim_end = curr_node.start + curr_node.size;
        if(curr_node.size <= max_leaf_size) {
            continue;
        }

        std::vector<Bucket> buckets(num_buckets);
        int min_axis = 0;
        int min_split_id = 0;
        BBox min_bbox_half_1;
        BBox min_bbox_half_2;
        size_t min_half_1_num_prim = 0;
        size_t min_half_2_num_prim = 0;
        float min_SAH = FLT_MAX;
        
        BBox curr_node_bbox = curr_node.bbox;
        for(int axis = 0; axis < 3; axis++) {
            // reset the buckets
            for(int i = 0; i < num_buckets; i++) {
                buckets[i].bbox.reset();
                buckets[i].num_prim = 0;
                buckets[i].total_area = 0.0f;
            }
            
            // assign each primitive of the node to a bucket
            for(size_t i = prim_start; i < prim_end; i++) {
                Primitive& curr_prim = primitives[i];
                Bucket& B = buckets[compute_bucket(curr_node_bbox, curr_prim.bbox().center(), axis, num_buckets)];
                B.bbox.enclose(curr_prim.bbox());
                B.num_prim++;
            }

            // evaluate the SAH for each possible split
            for(size_t i = 1; i < num_buckets; i++) {
                BBox bbox_half_1 = BBox();
                BBox bbox_half_2 = BBox();
                size_t half_1_num_prim = 0;
                size_t half_2_num_prim = 0;

                for(size_t j = 0; j < num_buckets; j++) {
                    if (j < i) {
                        // all the buckets before i
                        bbox_half_1.enclose(buckets[j].bbox);
                        half_1_num_prim += buckets[j].num_prim;
                    } else {
                        // all the buckets after i
                        bbox_half_2.enclose(buckets[j].bbox);
                        half_2_num_prim += buckets[j].num_prim;
                    }
                }

                float curr_SAH = (float)half_1_num_prim * bbox_half_1.surface_area() + (float)half_2_num_prim * bbox_half_2.surface_area();

                // update min cost
                if(curr_SAH < min_SAH) {
                    min_SAH = curr_SAH;
                    min_axis = axis;
                    min_split_id = (int)i;
                    min_bbox_half_1 = bbox_half_1;
                    min_bbox_half_2 = bbox_half_2;
                    min_half_2_num_prim = half_2_num_prim;
                    min_half_1_num_prim = half_1_num_prim;
                }
            }
        }

        // reordering the primitives based on the split
        std::partition(primitives.begin() + prim_start, primitives.begin() + prim_end, 
        [curr_node_bbox, min_axis, min_split_id](Primitive& prim) {
            float bbox_min = curr_node_bbox.min[min_axis];
            float bbox_max = curr_node_bbox.max[min_axis];

            float pos = (prim.bbox().center()[min_axis] - bbox_min) / (bbox_max - bbox_min);
            int idx = (size_t)(pos * (float)num_buckets);
            return idx < min_split_id;
        });

        size_t left_child_id = new_node(min_bbox_half_1, prim_start, min_half_1_num_prim, 0, 0);
        size_t right_child_id = new_node(min_bbox_half_2, prim_start + min_half_1_num_prim, min_half_2_num_prim, 0, 0);
        bbox_stack.push(left_child_id);
        bbox_stack.push(right_child_id);
        nodes[curr_id].l = left_child_id;
        nodes[curr_id].r = right_child_id;
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
    // TODO (PathTracer): Task 3
    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive
    // in the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.
    Trace ret;
    Vec2 range;
    range.x = 0.0f;
    range.y = FLT_MAX;
    if(nodes[root_idx].bbox.hit(ray, range)) {
        // iterate until leaf
        ret = find_closest_hit(ray, root_idx, range);
    }
    return ret;
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
