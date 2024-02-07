
#include <queue>
#include <set>
#include <unordered_map>
#include <iostream>

#include "../geometry/halfedge.h"
#include "debug.h"

/* Note on local operation return types:

    The local operations all return a std::optional<T> type. This is used so that your
    implementation can signify that it does not want to perform the operation for
    whatever reason (e.g. you don't want to allow the user to erase the last vertex).

    An optional can have two values: std::nullopt, or a value of the type it is
    parameterized on. In this way, it's similar to a pointer, but has two advantages:
    the value it holds need not be allocated elsewhere, and it provides an API that
    forces the user to check if it is null before using the value.

    In your implementaiton, if you have successfully performed the operation, you can
    simply return the required reference:

            ... collapse the edge ...
            return collapsed_vertex_ref;

    And if you wish to deny the operation, you can return the null optional:

            return std::nullopt;

    Note that the stubs below all reject their duties by returning the null optional.
*/

/*
    This method should replace the given vertex and all its neighboring
    edges and faces with a single face, returning the new face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::erase_vertex(Halfedge_Mesh::VertexRef v) {

    if(v->on_boundary()) {
        (void)v;
        return std::nullopt;
    }
    Halfedge_Mesh::FaceRef f = new_face();
    Halfedge_Mesh::HalfedgeRef h = v->halfedge();
    std::vector<Halfedge_Mesh::HalfedgeRef> h_twins;

    // erase all the halfedges and edges connected to the vertex
    do {
        erase(h);
        h_twins.push_back(h->twin());
        erase(h->twin());
        erase(h->edge());
        erase(h->face());
        h = h->twin()->next();
    } while(h != v->halfedge());
    erase(v);

    // walk along the contour of the new face
    std::reverse(h_twins.begin(), h_twins.end());
    h = v->halfedge()->next();
    int counter = 0;
    do {
        Halfedge_Mesh::HalfedgeRef h_before_v = h_twins[counter % h_twins.size()];

        while(h->next() != h_before_v) {
            h->face() = f;
            h = h->next();
        }

        Halfedge_Mesh::HalfedgeRef next_h = h_before_v->twin()->next();
        next_h->vertex()->halfedge() = next_h;
        h->next() = next_h;
        h->face() = f;
        h = next_h;
        counter++;

    } while(h != v->halfedge()->next());

    f->halfedge() = h;

    return f;
}


// helper functions
Halfedge_Mesh::HalfedgeRef find_prev_h(Halfedge_Mesh::HalfedgeRef h) {
    Halfedge_Mesh::HalfedgeRef h_prev = h;
    while(h_prev->next() != h) {
        h_prev = h_prev->next();
    }
    return h_prev;
}

/*
    This method should erase the given edge and return an iterator to the
    merged face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::erase_edge(Halfedge_Mesh::EdgeRef e) {

    if(e->on_boundary()) {
        (void) e;
        return std::nullopt;
    }

    Halfedge_Mesh::HalfedgeRef curr_h = e->halfedge();
    Halfedge_Mesh::FaceRef curr_f = curr_h->face(); // we keeping curr_f
    Halfedge_Mesh::HalfedgeRef curr_h_twin = curr_h->twin();
    Halfedge_Mesh::FaceRef oppo_f = curr_h_twin->face();

    Halfedge_Mesh::HalfedgeRef curr_h_prev = find_prev_h(curr_h);
    Halfedge_Mesh::HalfedgeRef curr_h_twin_prev = find_prev_h(curr_h_twin);

    Halfedge_Mesh::HalfedgeRef curr_h_next = curr_h->next();
    Halfedge_Mesh::HalfedgeRef curr_h_twin_next = curr_h_twin->next();

    if(curr_f == oppo_f) {
        (void) e;
        return std::nullopt;
    }

    Halfedge_Mesh::VertexRef curr_h_v = curr_h->vertex();
    Halfedge_Mesh::VertexRef curr_h_twin_v = curr_h_twin->vertex();

    // modify the halfedge of two vertices of deleted edge
    curr_h_v->halfedge() = curr_h_twin->next();
    curr_h_twin_v->halfedge() = curr_h->next();

    // modify face0's halfedge to something not deleted.
    curr_f->halfedge() = curr_h_twin->next();

    // loop through the halfedges of the opposite face and assign them to the other face
    Halfedge_Mesh::HalfedgeRef temp_h = curr_h_twin_next;
    while(temp_h != curr_h_twin) {
        temp_h->face() = curr_f;
        temp_h = temp_h->next();
    }

    // point curr_h_twin->prev and assign it's next to curr_h->next
    curr_h_twin_prev->next() = curr_h_next;

    // point curr_h->prev and assign it's next to curr_h_twin->next
    curr_h_prev->next() = curr_h_twin_next;

    erase(e);
    erase(oppo_f);
    erase(curr_h);
    erase(curr_h_twin);
    return curr_f;
}

bool is_rest_on_boundary(Halfedge_Mesh::HalfedgeRef h) {
    return h->next()->edge()->on_boundary() || h->next()->next()->edge()->on_boundary();
}

void collapse_non_triangle_face(Halfedge_Mesh::HalfedgeRef h, Halfedge_Mesh::VertexRef v) {
    Halfedge_Mesh::HalfedgeRef h_1 = h->next();
    Halfedge_Mesh::HalfedgeRef h_2 = find_prev_h(h);

    h_2->next() = h_1;
    h_1->vertex() = v;
    h_1->face()->halfedge() = h_1;
}

void collapse_triangle_face(Halfedge_Mesh::HalfedgeRef h, bool is_left, Halfedge_Mesh* mesh) {
    // the left traingle will disppear after the collapse
    Halfedge_Mesh::EdgeRef e_1 = h->next()->edge();
    Halfedge_Mesh::EdgeRef e_2 = h->next()->next()->edge();

    Halfedge_Mesh::HalfedgeRef h_1 = h->next();
    Halfedge_Mesh::HalfedgeRef h_2 = h_1->next();

    Halfedge_Mesh::HalfedgeRef h_1_twin = h_1->twin();
    Halfedge_Mesh::HalfedgeRef h_2_twin = h_2->twin();
    
    // these twins will become twins of each other
    h_1_twin->twin() = h_2_twin;
    h_2_twin->twin() = h_1_twin;

    Halfedge_Mesh::VertexRef v_1 = h_2->vertex();
    v_1->halfedge() = h_2_twin->next();

    if(is_left) {
        h_1_twin->edge() = e_2;
        e_2->halfedge() = h_2_twin;
        mesh->erase(e_1);
    } else {
        h_2_twin->edge() = e_1;
        e_1->halfedge() = h_1_twin;
        mesh->erase(e_2);
    }

    mesh->erase(h_1);
    mesh->erase(h_2);
}

/*
    This method should collapse the given edge and return an iterator to
    the new vertex created by the collapse.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(Halfedge_Mesh::EdgeRef e) {

    // if edge is boundary, return its vertex
    if (e->on_boundary()) {
        return e->halfedge()->vertex();
    }
    
    Halfedge_Mesh::HalfedgeRef curr_h = e->halfedge();
    Halfedge_Mesh::HalfedgeRef oppo_h = curr_h->twin();
    Halfedge_Mesh::HalfedgeRef input_h = e->halfedge();

    Halfedge_Mesh::FaceRef curr_f = curr_h->face();
    Halfedge_Mesh::FaceRef oppo_f = oppo_h->face();

    // check if either side of the edge is a triangle
    bool is_left_triangle = curr_f->degree() == 3;
    bool is_right_triangle = oppo_f->degree() == 3;

    // set curr vertex to the center of the edge
    Halfedge_Mesh::VertexRef curr_v = curr_h->vertex();
    Halfedge_Mesh::VertexRef oppo_v = oppo_h->vertex();
    curr_v->pos = e->center(); // we are keeping this vertex
    
    // update the vertex to point to the next next halfedge
    curr_v->halfedge() = curr_h->twin()->next()->twin()->next();
    oppo_v->halfedge() = oppo_h->twin()->next()->twin()->next();

    // update the left side
    if (is_left_triangle) {
        if (is_rest_on_boundary(curr_h)) {
            return e->halfedge()->vertex();
        }

        collapse_triangle_face(curr_h, true, this);

    } else {
        collapse_non_triangle_face(curr_h, curr_v);
    }

    // Since we are keeping the current vertex, we need to update the verticies of the edges that are connected to the top of the current edge
    curr_h = curr_h->next()->twin()->next();
    Halfedge_Mesh::HalfedgeRef h_to_update = curr_h;
    while(h_to_update != oppo_h) {
        h_to_update->vertex() = curr_v;
        h_to_update = h_to_update->twin()->next();
    }

    // update the right side
    if (is_right_triangle) {
        if(is_rest_on_boundary(oppo_h)) {
            return e->halfedge()->vertex();
        }

        collapse_triangle_face(oppo_h, false, this);

    } else {
        collapse_non_triangle_face(oppo_h, curr_v);
    }

    erase(input_h);
    erase(oppo_h);
    erase(e);
    erase(oppo_v);

    if(is_left_triangle) {
        erase(curr_f);
    }

    if(is_right_triangle) {
        erase(oppo_f);
    }

    return curr_v;
}

/*
    This method should collapse the given face and return an iterator to
    the new vertex created by the collapse.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_face(Halfedge_Mesh::FaceRef f) {

    (void)f;
    return std::nullopt;
}

/*
    This method should flip the given edge and return an iterator to the
    flipped edge.
*/
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::flip_edge(Halfedge_Mesh::EdgeRef e) {

    // Iterate half edge face and find the vertex in the middle of the edge path.
    HalfedgeRef he1 = e->halfedge();
    HalfedgeRef he1_fast = e->halfedge();

    HalfedgeRef curr = e->halfedge();
    do {
        curr = curr->next();
    } while (curr->id() != e->halfedge()->id());
    curr = e->halfedge()->twin();
    do {
        curr = curr->next();
    } while (curr->id() != e->halfedge()->twin()->id());

    HalfedgeRef new_he = new_halfedge();
    HalfedgeRef new_he_twin = new_halfedge();
    new_he->_twin = new_he_twin;
    new_he_twin->_twin = new_he;
    EdgeRef new_e = new_edge();
    new_e->_halfedge = new_he;
    new_he->_edge = new_e;
    new_he_twin->_edge = new_e;
    new_he->_face = he1->face();
    new_he_twin->_face = he1->twin()->face();

    // While the faster half edge reference has not completed one traversal around
    // the entire face, keep iterating.
    do {
        he1 = he1->next();
        he1_fast = he1_fast->next()->next();
    } while (he1_fast->id() != e->halfedge()->id() && he1_fast->next()->id() != e->halfedge()->id());
    HalfedgeRef he2 = he1->next();

    // Repeat for other face.
    HalfedgeRef he3 = e->halfedge()->twin();
    HalfedgeRef he3_fast = e->halfedge()->twin();
    do {
        he3 = he3->next();
        he3_fast = he3_fast->next()->next();
    } while (he3_fast->id() != e->halfedge()->twin()->id() && he3_fast->next()->id() != e->halfedge()->twin()->id());
    HalfedgeRef he4 = he3->next();

    // Also get the halfedges pointing to the edge we want to erase later.
    HalfedgeRef he5 = e->halfedge();
    while (he5->next()->id() != e->halfedge()->id()) {
        he5 = he5->next();
    }
    HalfedgeRef he6 = e->halfedge()->twin();
     while (he6->next()->id() != e->halfedge()->twin()->id()) {
        he6 = he6->next();
    }

    // The half-edge we want to return should point from the he2 vertex to the he4 vertex.
    // Update pointers.
    he1->next() = new_he;
    new_he->next() = he4;
    new_he->vertex() = he2->vertex();
    new_he->vertex()->_halfedge = new_he;

    he3->next() = new_he_twin;
    new_he_twin->next() = he2;
    new_he_twin->vertex() = he4->vertex();
    new_he_twin->vertex()->_halfedge = new_he_twin;

    he5->next() = e->halfedge()->twin()->next();
    he6->next() = e->halfedge()->next();
    e->halfedge()->vertex()->_halfedge = he5->next();
    e->halfedge()->twin()->vertex()->_halfedge = he6->next();
        
    // Unify halfedges->face and set face->halfedge.
    HalfedgeRef a = new_he->next();
    do {
        a->_face = new_he->face();
        a = a->next();
    } while (a->id() != new_he->id());
    a->face()->_halfedge = new_he;

    HalfedgeRef b = new_he_twin->next();
    do {
        b->_face = new_he_twin->face();
        b = b->next();
    } while (b->id() != new_he_twin->id());
    b->face()->_halfedge = new_he_twin;

    curr = new_e->halfedge();
    do {
        curr = curr->next();
    } while (curr->id() != new_e->halfedge()->id());

    curr = new_e->halfedge()->twin();
    do {
        curr = curr->next();
    } while (curr->id() != new_e->halfedge()->twin()->id());

    // Remove the original edge.
    erase(e);
    erase(e->halfedge());
    erase(e->halfedge()->twin());
    // Return a pointer to the new flipped edge.
    return new_e;
}

/*
    This method should split the given edge and return an iterator to the
    newly inserted vertex. The halfedge of this vertex should point along
    the edge that was split, rather than the new edges.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::split_edge(Halfedge_Mesh::EdgeRef e) {

    // Store references for old halfedges.
    HalfedgeRef he_north = e->halfedge();
    HalfedgeRef he_northwest = he_north->next();
    HalfedgeRef he_southwest = he_northwest->next();
    HalfedgeRef he_southeast = he_north->twin()->next();
    HalfedgeRef he_northeast = he_southeast->next();

    // Create new vertex at center.
    VertexRef new_v = new_vertex();

    // Create three new edges and their halfedges (we will reuse the original 
    // edge being split and its halfedges for the remaining edge and halfedges).
    EdgeRef e_east = new_edge();
    EdgeRef e_south = new_edge();
    EdgeRef e_west = new_edge();
    HalfedgeRef he_east = new_halfedge();
    HalfedgeRef he_east_twin = new_halfedge();
    HalfedgeRef he_south = new_halfedge();
    HalfedgeRef he_south_twin = new_halfedge();
    HalfedgeRef he_west = new_halfedge();
    HalfedgeRef he_west_twin = new_halfedge();


    // Create two new faces (we will reuse the original two faces for the remaining
    // two faces).
    FaceRef f_southeast = new_face();
    FaceRef f_southwest=  new_face();

    // Attach halfedges for the new vertex.
    new_v->_halfedge = he_north;
    new_v->pos = e->center(); // New vertex goes at the center of the edge.
    // Fix halfedge for old vertex.
    he_north->vertex()->_halfedge = he_southeast;

    // Attach halfedges for new and old faces.
    f_southeast->_halfedge = he_southeast;
    f_southwest->_halfedge = he_southwest;
    he_north->face()->_halfedge = he_north;
    he_north->twin()->face()->_halfedge = he_north->twin();
    he_southeast->_face = f_southeast;
    he_southwest->_face = f_southwest;

    // Attach halfedges to new edges.
    e_east->_halfedge = he_east;
    e_south->_halfedge = he_south;
    e_west->_halfedge = he_west;

    // Assign new halfedge information.
    he_north->_next = he_northwest;
    he_north->_vertex = new_v;
    he_north->twin()->_next = he_east;

    he_northeast->_next = he_north->twin();
    he_southeast->_next = he_east_twin;
    he_southwest->_next = he_south_twin;
    he_northwest->_next = he_west_twin;

    he_east->_next = he_northeast;
    he_east->_twin = he_east_twin;
    he_east->_vertex = new_v;
    he_east->_edge = e_east;
    he_east->_face = he_north->twin()->face();

    he_east_twin->_next = he_south;
    he_east_twin->_twin = he_east;
    he_east_twin->_vertex = he_northeast->vertex();
    he_east_twin->_edge = e_east;
    he_east_twin->_face = f_southeast;

    he_south->_next = he_southeast;
    he_south->_twin = he_south_twin;
    he_south->_vertex = new_v;
    he_south->_edge = e_south;
    he_south->_face = f_southeast;

    he_south_twin->_next = he_west;
    he_south_twin->_twin = he_south;
    he_south_twin->_vertex = he_southeast->vertex();
    he_south_twin->_edge = e_south;
    he_south_twin->_face = f_southwest;

    he_west->_next = he_southwest;
    he_west->_twin = he_west_twin;
    he_west->_vertex = new_v;
    he_west->_edge = e_west;
    he_west->_face = f_southwest;

    he_west_twin->_next = he_north;
    he_west_twin->_twin = he_west;
    he_west_twin->_vertex = he_southwest->vertex();
    he_west_twin->_edge = e_west;
    he_west_twin->_face = he_north->face();

    return new_v;
}

/* Note on the beveling process:

    Each of the bevel_vertex, bevel_edge, and bevel_face functions do not represent
    a full bevel operation. Instead, they should update the _connectivity_ of
    the mesh, _not_ the positions of newly created vertices. In fact, you should set
    the positions of new vertices to be exactly the same as wherever they "started from."

    When you click on a mesh element while in bevel mode, one of those three functions
    is called. But, because you may then adjust the distance/offset of the newly
    beveled face, we need another method of updating the positions of the new vertices.

    This is where bevel_vertex_positions, bevel_edge_positions, and
    bevel_face_positions come in: these functions are called repeatedly as you
    move your mouse, the position of which determins the normal and tangent offset
    parameters. These functions are also passed an array of the original vertex
    positions: for  bevel_vertex, it has one element, the original vertex position,
    for bevel_edge,  two for the two vertices, and for bevel_face, it has the original
    position of each vertex in halfedge order. You should use these positions, as well
    as the normal and tangent offset fields to assign positions to the new vertices.

    Finally, note that the normal and tangent offsets are not relative values - you
    should compute a particular new position from them, not a delta to apply.
*/

/*
    This method should replace the vertex v with a face, corresponding to
    a bevel operation. It should return the new face.  NOTE: This method is
    responsible for updating the *connectivity* of the mesh only---it does not
    need to update the vertex positions.  These positions will be updated in
    Halfedge_Mesh::bevel_vertex_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_vertex(Halfedge_Mesh::VertexRef v) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."

    (void)v;
    return std::nullopt;
}

/*
    This method should replace the edge e with a face, corresponding to a
    bevel operation. It should return the new face. NOTE: This method is
    responsible for updating the *connectivity* of the mesh only---it does not
    need to update the vertex positions.  These positions will be updated in
    Halfedge_Mesh::bevel_edge_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_edge(Halfedge_Mesh::EdgeRef e) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."

    (void)e;
    return std::nullopt;
}

/*
    This method should replace the face f with an additional, inset face
    (and ring of faces around it), corresponding to a bevel operation. It
    should return the new face.  NOTE: This method is responsible for updating
    the *connectivity* of the mesh only---it does not need to update the vertex
    positions. These positions will be updated in
    Halfedge_Mesh::bevel_face_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_face(Halfedge_Mesh::FaceRef f) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."

    // edge case - can't bevel boundary face
    if(f->is_boundary()) {
        (void) f;
        return std::nullopt;
    }

    Halfedge_Mesh::HalfedgeRef h = f->halfedge();
    Halfedge_Mesh::HalfedgeRef h_placeholder = h;
    int degree = f->degree();

    std::vector<Halfedge_Mesh::VertexRef> new_vertices;
    std::vector<Halfedge_Mesh::EdgeRef> new_edges;
    std::vector<Halfedge_Mesh::FaceRef> new_faces;
    std::vector<Halfedge_Mesh::HalfedgeRef> new_halfedges;

    // allocate new components. we create a new face, a new vertex, two new edges
    // and four new edges for every side of the input face
    // the two new edges are the left and bottom edges of the new face
    for (int i = 0; i < degree; i++) {
        new_vertices.push_back(new_vertex());
        new_faces.push_back(new_face());
    }

    for (int i = 0; i < degree * 2; i++) {
        new_edges.push_back(new_edge());
    }

    for (int i = 0; i < degree * 4; i++) {
        new_halfedges.push_back(new_halfedge());
    }

    f->halfedge() = new_halfedges[3];

    for(int i = 0; i < degree; ++i) {
        h = h_placeholder;

        Halfedge_Mesh::VertexRef curr_v = new_vertices[i];
        Halfedge_Mesh::FaceRef curr_f = new_faces[i];
        Halfedge_Mesh::EdgeRef curr_e0 = new_edges[i * 2];
        Halfedge_Mesh::EdgeRef curr_e1 = new_edges[i * 2 + 1];
        // h0 is the next halfedge of h
        Halfedge_Mesh::HalfedgeRef curr_h0 = new_halfedges[i * 4];
        // h1 is the next halfedge of h0
        Halfedge_Mesh::HalfedgeRef curr_h1 = new_halfedges[i * 4 + 1];
        // h2 is the next halfedge of h1 
        Halfedge_Mesh::HalfedgeRef curr_h2 = new_halfedges[i * 4 + 2];
        // h3 is the twin of h1
        Halfedge_Mesh::HalfedgeRef curr_h3 = new_halfedges[i * 4 + 3];
        
        // top left coner of the new face
        Halfedge_Mesh::VertexRef top_left_corner = h->next()->vertex();

        // h now points to the new face
        h->face() = curr_f;

        // new vertex is the bottom left corner of the new face
        curr_v->halfedge() = curr_h1;
        curr_v->pos = top_left_corner->pos; // pos is the vertex that h points to

        // new edges' halfedges
        curr_e0->halfedge() = curr_h0;
        curr_e1->halfedge() = curr_h1;

        // new halfedges' halfedges
        curr_f->halfedge() = curr_h1;

        // half edges
        int next_h2_index = (i + 2) * 4 - 2;
        int num_h = degree * 4;
        Halfedge_Mesh::HalfedgeRef next_h2 = new_halfedges[next_h2_index % num_h];

        curr_h0->Halfedge_Mesh::Halfedge::set_neighbors(curr_h1, next_h2, top_left_corner, curr_e0, curr_f);


        // h1 follows h0
        curr_h1->Halfedge_Mesh::Halfedge::set_neighbors(curr_h2, curr_h3, curr_v, curr_e1, curr_f);


        // h2 follows h1
        int prev_h0_index = (i + degree - 1) * 4; // twin of h2
        Halfedge_Mesh::HalfedgeRef prev_h0 = new_halfedges[prev_h0_index % num_h];

        int prev_edge_index = (i + degree - 1) * 2; // edge of h2
        Halfedge_Mesh::EdgeRef prev_edge = new_edges[prev_edge_index % (degree * 2)];

        int prev_vertex_index = (i + degree - 1) % degree;
        Halfedge_Mesh::VertexRef prev_vertex = new_vertices[prev_vertex_index];

        curr_h2->Halfedge_Mesh::Halfedge::set_neighbors(h, prev_h0, prev_vertex, prev_edge, curr_f);


        // h3 is the twin of h1, all the h3s form the halfedges of the beveled face
        int next_h3_index = (i + 2) * 4 - 1; // the next of h3
        Halfedge_Mesh::HalfedgeRef next_h3 = new_halfedges[next_h3_index % num_h];

        curr_h3->Halfedge_Mesh::Halfedge::set_neighbors(next_h3, curr_h1, prev_vertex, curr_e1, f);

       
        h_placeholder = h_placeholder->next();
        h->next() = curr_h0;
    }
    
    return f;
}



/*
    Compute new vertex positions for the vertices of the beveled vertex.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the original vertex position and its associated outgoing edge
    to compute a new vertex position along the outgoing edge.
*/
void Halfedge_Mesh::bevel_vertex_positions(const std::vector<Vec3>& start_positions,
                                           Halfedge_Mesh::FaceRef face, float tangent_offset) {

    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());

    (void)new_halfedges;
    (void)start_positions;
    (void)face;
    (void)tangent_offset;
}

/*
    Compute new vertex positions for the vertices of the beveled edge.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the preceding and next vertex position from the original mesh
    (in the orig array) to compute an offset vertex position.

    Note that there is a 1-to-1 correspondence between halfedges in
    newHalfedges and vertex positions
    in orig.  So, you can write loops of the form

    for(size_t i = 0; i < new_halfedges.size(); i++)
    {
            Vector3D pi = start_positions[i]; // get the original vertex
            position corresponding to vertex i
    }
*/
void Halfedge_Mesh::bevel_edge_positions(const std::vector<Vec3>& start_positions,
                                         Halfedge_Mesh::FaceRef face, float tangent_offset) {

    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());

    (void)new_halfedges;
    (void)start_positions;
    (void)face;
    (void)tangent_offset;
}

/*
    Compute new vertex positions for the vertices of the beveled face.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the preceding and next vertex position from the original mesh
    (in the start_positions array) to compute an offset vertex
    position.

    Note that there is a 1-to-1 correspondence between halfedges in
    new_halfedges and vertex positions
    in orig. So, you can write loops of the form

    for(size_t i = 0; i < new_halfedges.size(); i++)
    {
            Vec3 pi = start_positions[i]; // get the original vertex
            position corresponding to vertex i
    }
*/
void Halfedge_Mesh::bevel_face_positions(const std::vector<Vec3>& start_positions,
                                         Halfedge_Mesh::FaceRef face, float tangent_offset,
                                         float normal_offset) {

    if(flip_orientation) normal_offset = -normal_offset;
    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());

    (void)new_halfedges;
    (void)start_positions;
    (void)face;
    (void)tangent_offset;
    (void)normal_offset;

    int new_halfedges_size = new_halfedges.size();
    Vec3 face_normal = face->normal();

    for(int i = 0; i < new_halfedges_size; i++) {
        Vec3 pi = start_positions[i];
        Vec3 pi_prev = start_positions[(i + new_halfedges_size - 1) % new_halfedges_size];
        Vec3 pi_next = start_positions[(i + 1) % new_halfedges_size];

        Vec3 normal_offset_vector = face_normal * normal_offset;
        Vec3 tangent_offset_vector = (pi - (pi_prev + pi_next) / 2.f).unit() * tangent_offset;

        Vec3 new_pos = pi + tangent_offset_vector - normal_offset_vector;
        new_halfedges[i]->vertex()->pos = new_pos;
    }
}

/*
    Splits all non-triangular faces into triangles.
*/
void Halfedge_Mesh::triangulate() {
    // For each face, choose a vertex u. Iterate through all other vertices v around
    // the face and add an edge u -> v if there doesn't already exist one.
    for (FaceRef curr_face = faces_begin(); curr_face != faces_end(); curr_face++) {
        // Store all half edges and vertices around this face.
        std::vector<HalfedgeRef> face_hes;
        std::vector<VertexRef> vs;
        HalfedgeRef he = curr_face->halfedge();
        do {
            face_hes.push_back(he);
            vs.push_back(he->vertex());
            he = he->next();
        } while (he != curr_face->halfedge());
        if ((int)vs.size() <= 3) continue; // Don't need to triangulate a triangle.
        VertexRef base_v = vs[0]; // This is the vertex out of which all new edges should come.

        // Create all new faces, edges and halfedges.
        std::vector<FaceRef> new_faces;
        std::vector<EdgeRef> new_edges;
        std::vector<HalfedgeRef> new_halfedges;
        std::vector<HalfedgeRef> new_halfedge_twins;
        for (int i = 0; i < (int)face_hes.size() - 3; i++) {
            HalfedgeRef new_he = new_halfedge();
            HalfedgeRef new_he_twin = new_halfedge();
            EdgeRef new_e = new_edge();
            FaceRef new_f = new_face();

            // Attach all halfedges to edges/vertices, edges to halfedge, and faces to halfedges.
            new_he->_twin = new_he_twin;
            new_he->_vertex = base_v;
            new_he->_edge = new_e;
            new_he_twin->_twin = new_he;
            new_he_twin->_vertex = vs[i + 2];
            new_he_twin->_edge = new_e;
            new_e->_halfedge = new_he;
            new_f->_halfedge = new_he_twin;
            
            new_faces.push_back(new_f);
            new_edges.push_back(new_e);
            new_halfedges.push_back(new_he);
            new_halfedge_twins.push_back(new_he_twin);
        }

        // Attach each new halfedge to its face.
        for (int i = 0; i < (int)new_faces.size(); i++) {
            new_halfedge_twins[i]->_face = new_faces[i];
            if (i < (int)new_faces.size() - 1) {
                new_halfedges[i]->_face = new_faces[i + 1];
            } else {
                new_halfedges[i]->_face = curr_face;
            }
        }
        // Reassign the halfedge of the original face.
        // face_hes[face_hes.size() - 1]->_face = curr_face;

        // Reassign faces of original face halfedges.
        for (int i = 0; i < (int)face_hes.size(); i++) {
            if (i < 1) {
                face_hes[i]->_face = new_faces[0];
            } else if (i > (int)face_hes.size() - 3) {
                face_hes[i]->_face = curr_face;
            } else {
                face_hes[i]->_face = new_faces[i - 1];
            }
        }
        curr_face->_halfedge = face_hes[face_hes.size() - 1];

        // Assign nexts to new halfedges.
        for (int i = 0; i < (int)new_halfedges.size(); i++) {
            new_halfedges[i]->_next = face_hes[i + 2];
            if (i > 0) {
                new_halfedge_twins[i]->_next = new_halfedges[i - 1];
            } else {
                new_halfedge_twins[i]->_next = face_hes[0];
            }
            face_hes[i + 1]->_next = new_halfedge_twins[i];
        }
        // Reassgn next of last old face halfedge.
        face_hes[face_hes.size() - 1]->_next = new_halfedges[new_halfedges.size() - 1];
    }
}

/* Note on the quad subdivision process:

        Unlike the local mesh operations (like bevel or edge flip), we will perform
        subdivision by splitting *all* faces into quads "simultaneously."  Rather
        than operating directly on the halfedge data structure (which as you've
        seen is quite difficult to maintain!) we are going to do something a bit nicer:
           1. Create a raw list of vertex positions and faces (rather than a full-
              blown halfedge mesh).
           2. Build a new halfedge mesh from these lists, replacing the old one.
        Sometimes rebuilding a data structure from scratch is simpler (and even
        more efficient) than incrementally modifying the existing one.  These steps are
        detailed below.

  Step I: Compute the vertex positions for the subdivided mesh.
        Here we're going to do something a little bit strange: since we will
        have one vertex in the subdivided mesh for each vertex, edge, and face in
        the original mesh, we can nicely store the new vertex *positions* as
        attributes on vertices, edges, and faces of the original mesh. These positions
        can then be conveniently copied into the new, subdivided mesh.
        This is what you will implement in linear_subdivide_positions() and
        catmullclark_subdivide_positions().

  Steps II-IV are provided (see Halfedge_Mesh::subdivide()), but are still detailed
  here:

  Step II: Assign a unique index (starting at 0) to each vertex, edge, and
        face in the original mesh. These indices will be the indices of the
        vertices in the new (subdivided mesh).  They do not have to be assigned
        in any particular order, so long as no index is shared by more than one
        mesh element, and the total number of indices is equal to V+E+F, i.e.,
        the total number of vertices plus edges plus faces in the original mesh.
        Basically we just need a one-to-one mapping between original mesh elements
        and subdivided mesh vertices.

  Step III: Build a list of quads in the new (subdivided) mesh, as tuples of
        the element indices defined above. In other words, each new quad should be
        of the form (i,j,k,l), where i,j,k and l are four of the indices stored on
        our original mesh elements.  Note that it is essential to get the orientation
        right here: (i,j,k,l) is not the same as (l,k,j,i).  Indices of new faces
        should circulate in the same direction as old faces (think about the right-hand
        rule).

  Step IV: Pass the list of vertices and quads to a routine that clears
        the internal data for this halfedge mesh, and builds new halfedge data from
        scratch, using the two lists.
*/

/*
    Compute new vertex positions for a mesh that splits each polygon
    into quads (by inserting a vertex at the face midpoint and each
    of the edge midpoints).  The new vertex positions will be stored
    in the members Vertex::new_pos, Edge::new_pos, and
    Face::new_pos.  The values of the positions are based on
    simple linear interpolation, e.g., the edge midpoints and face
    centroids.
*/
void Halfedge_Mesh::linear_subdivide_positions() {

    // For each vertex, assign Vertex::new_pos to
    // its original position, Vertex::pos.
    for (VertexRef v = vertices_begin(); v != vertices_end(); v++) {
        v->new_pos = v->pos;
    }

    // For each edge, assign the midpoint of the two original
    // positions to Edge::new_pos.
    for (EdgeRef e = edges_begin(); e != edges_end(); e++) {
        // VertexRef u = e->halfedge()->vertex();
        // VertexRef v = e->halfedge()->next()->vertex();
        // Vec3 midpoint = (u->pos + v->pos) / 2;
        e->new_pos = e->center();
    }

    // For each face, assign the centroid (i.e., arithmetic mean)
    // of the original vertex positions to Face::new_pos. Note
    // that in general, NOT all faces will be triangles!
    for (FaceRef f = faces_begin(); f != faces_end(); f++) {
        std::vector<Vec3> v_positions;
        HalfedgeRef he = f->halfedge();
        do {
            v_positions.push_back(he->vertex()->pos);
            he = he->next();
        } while (he != f->halfedge());

        Vec3 centroid = v_positions[0];
        for (int i = 1; i < (int)v_positions.size(); i++) {
            centroid += v_positions[i];
        }
        f->new_pos = centroid / (int)v_positions.size();
    }
}

/*
    Compute new vertex positions for a mesh that splits each polygon
    into quads (by inserting a vertex at the face midpoint and each
    of the edge midpoints).  The new vertex positions will be stored
    in the members Vertex::new_pos, Edge::new_pos, and
    Face::new_pos.  The values of the positions are based on
    the Catmull-Clark rules for subdivision.

    Note: this will only be called on meshes without boundary
*/
void Halfedge_Mesh::catmullclark_subdivide_positions() {

    // The implementation for this routine should be
    // a lot like Halfedge_Mesh:linear_subdivide_positions:(),
    // except that the calculation of the positions themsevles is
    // slightly more involved, using the Catmull-Clark subdivision
    // rules. (These rules are outlined in the Developer Manual.)

    // Faces
    for (FaceRef f = faces_begin(); f != faces_end(); f++) {
        std::vector<Vec3> v_positions;
        HalfedgeRef he = f->halfedge();
        do {
            v_positions.push_back(he->vertex()->pos);
            he = he->next();
        } while (he != f->halfedge());

        Vec3 centroid = v_positions[0];
        for (int i = 1; i < (int)v_positions.size(); i++) {
            centroid += v_positions[i];
        }
        f->new_pos = centroid / (int)v_positions.size();
    }

    // Edges
    for (EdgeRef e = edges_begin(); e != edges_end(); e++) {
        Vec3 new_edge_pos = Vec3(0, 0, 0);
        new_edge_pos += e->halfedge()->face()->new_pos;
        new_edge_pos += e->halfedge()->twin()->face()->new_pos;
        new_edge_pos += e->halfedge()->vertex()->pos;
        new_edge_pos += e->halfedge()->twin()->vertex()->pos;
        e->new_pos = new_edge_pos / 4;
    }

    // Vertices
    for (VertexRef v = vertices_begin(); v != vertices_end(); v++) {
        // Get all adjacent faces.
        HalfedgeRef start_he = v->halfedge();
        HalfedgeRef curr_he = start_he;
        std::vector<Vec3> adj_faces_new_pos;
        std::vector<Vec3> all_edges_midpoints;
        do {
            // Iterate around adjacent faces to hop to the next face via twin edges.
            adj_faces_new_pos.push_back(curr_he->face()->new_pos);
            Vec3 adj_edge_midpoint = curr_he->edge()->center();
            all_edges_midpoints.push_back(adj_edge_midpoint);
            do {
                curr_he = curr_he->next();
            } while (curr_he->next() != start_he);
            start_he = curr_he->twin();
            curr_he = start_he;
        } while (curr_he != v->halfedge());

        Vec3 q = Vec3(0, 0, 0);
        for (int i = 0; i < (int)adj_faces_new_pos.size(); i++) {
            q += adj_faces_new_pos[i];
        }
        q /= (int)adj_faces_new_pos.size();

        Vec3 r = Vec3(0, 0, 0);
        for (int i = 0; i < (int)all_edges_midpoints.size(); i++) {
            r += all_edges_midpoints[i];
        }
        r /= (int)all_edges_midpoints.size();

        Vec3 s = v->pos;
        int n = (int)adj_faces_new_pos.size();
        v->new_pos = (q + 2 * r + (n - 3) * s) / n;
    }
}

/*
        This routine should increase the number of triangles in the mesh
        using Loop subdivision. Note: this is will only be called on triangle meshes.
*/
void Halfedge_Mesh::loop_subdivide() {

    // Compute new positions for all the vertices in the input mesh, using
    // the Loop subdivision rule, and store them in Vertex::new_pos.
    // -> At this point, we also want to mark each vertex as being a vertex of the
    //    original mesh. Use Vertex::is_new for this.
    // -> Next, compute the updated vertex positions associated with edges, and
    //    store it in Edge::new_pos.
    // -> Next, we're going to split every edge in the mesh, in any order.  For
    //    future reference, we're also going to store some information about which
    //    subdivided edges come from splitting an edge in the original mesh, and
    //    which edges are new, by setting the flat Edge::is_new. Note that in this
    //    loop, we only want to iterate over edges of the original mesh.
    //    Otherwise, we'll end up splitting edges that we just split (and the
    //    loop will never end!)
    // -> Now flip any new edge that connects an old and new vertex.
    // -> Finally, copy the new vertex positions into final Vertex::pos.

    // Each vertex and edge of the original surface can be associated with a
    // vertex in the new (subdivided) surface.
    // Therefore, our strategy for computing the subdivided vertex locations is to
    // *first* compute the new positions
    // using the connectivity of the original (coarse) mesh; navigating this mesh
    // will be much easier than navigating
    // the new subdivided (fine) mesh, which has more elements to traverse.  We
    // will then assign vertex positions in
    // the new mesh based on the values we computed for the original mesh.

    // Compute updated positions for all the vertices in the original mesh, using
    // the Loop subdivision rule.

    // Next, compute the updated vertex positions associated with edges.

    // Next, we're going to split every edge in the mesh, in any order. For
    // future reference, we're also going to store some information about which
    // subdivided edges come from splitting an edge in the original mesh, and
    // which edges are new.
    // In this loop, we only want to iterate over edges of the original
    // mesh---otherwise, we'll end up splitting edges that we just split (and
    // the loop will never end!)

    // Finally, flip any new edge that connects an old and new vertex.

    // Copy the updated vertex positions to the subdivided mesh.
}

/*
    Isotropic remeshing. Note that this function returns success in a similar
    manner to the local operations, except with only a boolean value.
    (e.g. you may want to return false if this is not a triangle mesh)
*/
bool Halfedge_Mesh::isotropic_remesh() {

    // Compute the mean edge length.
    // Repeat the four main steps for 5 or 6 iterations
    // -> Split edges much longer than the target length (being careful about
    //    how the loop is written!)
    // -> Collapse edges much shorter than the target length.  Here we need to
    //    be EXTRA careful about advancing the loop, because many edges may have
    //    been destroyed by a collapse (which ones?)
    // -> Now flip each edge if it improves vertex degree
    // -> Finally, apply some tangential smoothing to the vertex positions

    // Note: if you erase elements in a local operation, they will not be actually deleted
    // until do_erase or validate are called. This is to facilitate checking
    // for dangling references to elements that will be erased.
    // The rest of the codebase will automatically call validate() after each op,
    // but here simply calling collapse_edge() will not erase the elements.
    // You should use collapse_edge_erase() instead for the desired behavior.

    return false;
}

/* Helper type for quadric simplification */
struct Edge_Record {
    Edge_Record() {
    }
    Edge_Record(std::unordered_map<Halfedge_Mesh::VertexRef, Mat4>& vertex_quadrics,
                Halfedge_Mesh::EdgeRef e)
        : edge(e) {

        // Compute the combined quadric from the edge endpoints.
        Halfedge_Mesh::VertexRef u = e->halfedge()->vertex();
        Halfedge_Mesh::VertexRef v = e->halfedge()->twin()->vertex();
        // -> Build the 3x3 linear system whose solution minimizes the quadric error
        //    associated with these two endpoints.
        Mat4 K_sum = vertex_quadrics[u] + vertex_quadrics[v];
        Mat4 A(K_sum);
        Vec3 b = -K_sum[3].xyz();
        for(int i = 0; i < 3; ++i) {
            A[3][i] = A[i][3] = 0.0f;
        }
        A[3][3] = 1.0f;

        // -> Use this system to solve for the optimal position, and store it in
        //    Edge_Record::optimal.
        optimal = e->center(); // non invertible case

        if (std::abs(A.det()) > 1e-8) {
            optimal = A.inverse() * b;
        }

        // -> Also store the cost associated with collapsing this edge in
        //    Edge_Record::cost.
        cost = dot(Vec4(optimal, 1.0f), (K_sum * Vec4(optimal, 1.0f)));
    }

    Halfedge_Mesh::EdgeRef edge;
    Vec3 optimal;
    float cost;
};

/* Comparison operator for Edge_Records so std::set will properly order them */
bool operator<(const Edge_Record& r1, const Edge_Record& r2) {
    if(r1.cost != r2.cost) {
        return r1.cost < r2.cost;
    }
    Halfedge_Mesh::EdgeRef e1 = r1.edge;
    Halfedge_Mesh::EdgeRef e2 = r2.edge;
    return &*e1 < &*e2;
}

/** Helper type for quadric simplification
 *
 * A PQueue is a minimum-priority queue that
 * allows elements to be both inserted and removed from the
 * queue.  Together, one can easily change the priority of
 * an item by removing it, and re-inserting the same item
 * but with a different priority.  A priority queue, for
 * those who don't remember or haven't seen it before, is a
 * data structure that always keeps track of the item with
 * the smallest priority or "score," even as new elements
 * are inserted and removed.  Priority queues are often an
 * essential component of greedy algorithms, where one wants
 * to iteratively operate on the current "best" element.
 *
 * PQueue is templated on the type T of the object
 * being queued.  For this reason, T must define a comparison
 * operator of the form
 *
 *    bool operator<( const T& t1, const T& t2 )
 *
 * which returns true if and only if t1 is considered to have a
 * lower priority than t2.
 *
 * Basic use of a PQueue might look
 * something like this:
 *
 *    // initialize an empty queue
 *    PQueue<myItemType> queue;
 *
 *    // add some items (which we assume have been created
 *    // elsewhere, each of which has its priority stored as
 *    // some kind of internal member variable)
 *    queue.insert( item1 );
 *    queue.insert( item2 );
 *    queue.insert( item3 );
 *
 *    // get the highest priority item currently in the queue
 *    myItemType highestPriorityItem = queue.top();
 *
 *    // remove the highest priority item, automatically
 *    // promoting the next-highest priority item to the top
 *    queue.pop();
 *
 *    myItemType nextHighestPriorityItem = queue.top();
 *
 *    // Etc.
 *
 *    // We can also remove an item, making sure it is no
 *    // longer in the queue (note that this item may already
 *    // have been removed, if it was the 1st or 2nd-highest
 *    // priority item!)
 *    queue.remove( item2 );
 *
 */
template<class T> struct PQueue {
    void insert(const T& item) {
        queue.insert(item);
    }
    void remove(const T& item) {
        if(queue.find(item) != queue.end()) {
            queue.erase(item);
        }
    }
    const T& top(void) const {
        return *(queue.begin());
    }
    void pop(void) {
        queue.erase(queue.begin());
    }
    size_t size() {
        return queue.size();
    }

    std::set<T> queue;
};


void restore_edge_information(Halfedge_Mesh::VertexRef v,
                              std::unordered_map<Halfedge_Mesh::EdgeRef, Edge_Record>& edge_records,
                              PQueue<Edge_Record>& edge_queue,
                              std::unordered_map<Halfedge_Mesh::EdgeRef, bool>& edge_not_deleted,
                              std::unordered_map<Halfedge_Mesh::VertexRef, Mat4>& vertex_quadrics) {
    
    Halfedge_Mesh::HalfedgeRef h = v->halfedge();
    do {
        if(edge_not_deleted[h->edge()] != true) {
            Edge_Record new_record(vertex_quadrics, h->edge());
            edge_queue.insert(new_record);
            edge_records[h->edge()] = new_record;
        }
        h = h->twin()->next();
    } while(h != v->halfedge());
}

void remove_edge_information(Halfedge_Mesh::VertexRef v,
                             std::unordered_map<Halfedge_Mesh::EdgeRef, Edge_Record>& edge_records,
                             PQueue<Edge_Record>& edge_queue) {

    Halfedge_Mesh::HalfedgeRef h = v->halfedge();
    do {
        edge_queue.remove(edge_records[h->edge()]);
        h = h->twin()->next();
    } while(h != v->halfedge());
}

/*
    Mesh simplification. Note that this function returns success in a similar
    manner to the local operations, except with only a boolean value.
    (e.g. you may want to return false if you can't simplify the mesh any
    further without destroying it.)
*/
bool Halfedge_Mesh::simplify() {

    std::unordered_map<VertexRef, Mat4> vertex_quadrics;
    std::unordered_map<FaceRef, Mat4> face_quadrics;
    std::unordered_map<EdgeRef, Edge_Record> edge_records;
    PQueue<Edge_Record> edge_queue;

    // Compute initial quadrics for each face by simply writing the plane equation
    // for the face in homogeneous coordinates. These quadrics should be stored
    // in face_quadrics
    for(FaceRef f = faces_begin(); f != faces_end(); f++) {
        if(f->degree() != 3)
            return false;

        Vec3 n = f->normal();
        float d = -dot(n, f->halfedge()->vertex()->pos);
        Vec4 v(n, d);
        Mat4 K = outer(v, v);
        face_quadrics.insert({f, K});
    }

    // -> Compute an initial quadric for each vertex as the sum of the quadrics
    //    associated with the incident faces, storing it in vertex_quadrics
    for(VertexRef v = vertices_begin(); v != vertices_end(); v++) {
        Mat4 K(Mat4::Zero);
        HalfedgeRef h = v->halfedge();
        do {
            FaceRef curr_f = h->face();
            K += face_quadrics[curr_f];
            h = h->twin()->next();
        } while(h != v->halfedge());

        vertex_quadrics.insert({v, K});
    }


    // -> Build a priority queue of edges according to their quadric error cost,
    //    i.e., by building an Edge_Record for each edge and sticking it in the
    //    queue. You may want to use the above PQueue<Edge_Record> for this.
    std::unordered_map<EdgeRef, bool> edge_not_deleted;

    for(EdgeRef e = edges_begin(); e != edges_end(); e++) {
        Edge_Record curr_record(vertex_quadrics, e);
        edge_records.insert({e, curr_record});
        edge_queue.insert(curr_record);
        edge_not_deleted.insert({e, false});
    }

    // -> Until we reach the target edge budget, collapse the best edge. Remember
    //    to remove from the queue any edge that touches the collapsing edge
    //    BEFORE it gets collapsed, and add back into the queue any edge touching
    //    the collapsed vertex AFTER it's been collapsed. Also remember to assign
    //    a quadric to the collapsed vertex, and to pop the collapsed edge off the
    //    top of the queue.
    int input_face_count = (int)faces.size();
    int target_face_count = ceil(input_face_count / 4);
    std::cout << "input_face_count: " << input_face_count << std::endl;
    std::cout << "target_face_count: " << target_face_count << std::endl;

    while ((int)faces.size() > target_face_count && !edge_queue.queue.empty()) {
        std::cout << "face size: " << (int)faces.size() << std::endl;

        // 1. Get the cheapest edge from the queue.
        Edge_Record best_edge = edge_queue.top();

        // 2. Remove the cheapest edge from the queue by calling pop().
        edge_queue.pop();

        // 3. Compute the new quadric by summing the quadrics at its two endpoints.
        EdgeRef e = best_edge.edge;
        VertexRef v0 = e->halfedge()->vertex();
        VertexRef v1 = e->halfedge()->twin()->vertex();
        Mat4 new_K = vertex_quadrics[v0] + vertex_quadrics[v1];

        // 4. Remove any edge touching either of its endpoints from the queue.
        remove_edge_information(v0, edge_records, edge_queue);
        remove_edge_information(v1, edge_records, edge_queue);

        // 5. Collapse the edge.
        std::optional<Halfedge_Mesh::VertexRef> v_erased = collapse_edge_erase(e);

        if(v_erased.has_value() == false) {
            // if the edge is not deleted, we need to add the edge back to the queue
            edge_not_deleted[e] = true;
            restore_edge_information(v0, edge_records, edge_queue, edge_not_deleted, vertex_quadrics);

            restore_edge_information(v1, edge_records, edge_queue, edge_not_deleted, vertex_quadrics);
            continue;
        }

        // 6. Set the quadric of the new vertex to the quadric computed in Step 3.
        VertexRef new_v = v_erased.value();
        new_v->pos = best_edge.optimal;
        vertex_quadrics[new_v] = new_K;

        // 7. Insert any edge touching the new vertex into the queue, creating new edge records for each of them.
        Halfedge_Mesh::HalfedgeRef h = new_v->halfedge();
        do {
            Edge_Record new_record(vertex_quadrics, h->edge());
            edge_queue.insert(new_record);
            edge_records[h->edge()] = new_record;
            h = h->twin()->next();
        } while(h != new_v->halfedge());
    }


    // Note: if you erase elements in a local operation, they will not be actually deleted
    // until do_erase or validate are called. This is to facilitate checking
    // for dangling references to elements that will be erased.
    // The rest of the codebase will automatically call validate() after each op,
    // but here simply calling collapse_edge() will not erase the elements.
    // You should use collapse_edge_erase() instead for the desired behavior.

    return true;
}
