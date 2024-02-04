
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

    (void)v;
    return std::nullopt;
}

/*
    This method should erase the given edge and return an iterator to the
    merged face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::erase_edge(Halfedge_Mesh::EdgeRef e) {

    (void)e;
    return std::nullopt;
}

/*
    This method should collapse the given edge and return an iterator to
    the new vertex created by the collapse.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(Halfedge_Mesh::EdgeRef e) {

    (void)e;
    return std::nullopt;
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

    (void)e;
    // Iterate half edge face and find the vertex in the middle of the edge path.
    HalfedgeRef he1 = e->halfedge();
    HalfedgeRef he1_fast = e->halfedge();

    std::cout << "Original edge + face info:\n";
    HalfedgeRef curr = e->halfedge();
    do {
        std::cout << "he: " + std::to_string(curr->id()) + " v: " + std::to_string(curr->vertex()->id()) + " f: " + std::to_string(curr->face()->id()) + "\n";
        curr = curr->next();
    } while (curr->id() != e->halfedge()->id());
    std::cout << "Original twin edge + face info:\n";
    curr = e->halfedge()->twin();
    do {
        std::cout << "he: " + std::to_string(curr->id()) + " v: " + std::to_string(curr->vertex()->id()) + " f: " + std::to_string(curr->face()->id()) + "\n";
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

    std::cout << "New_e he: " + std::to_string(new_e->halfedge()->id()) + "\n";

    // While the faster half edge reference has not completed one traversal around
    // the entire face, keep iterating.
    do {
        he1 = he1->next();
        he1_fast = he1_fast->next()->next();
    } while (he1_fast->id() != e->halfedge()->id() && he1_fast->next()->id() != e->halfedge()->id());
    HalfedgeRef he2 = he1->next();
    std::cout << "he2 id: " + std::to_string(he2->id()) + "\n";

    // Repeat for other face.
    HalfedgeRef he3 = e->halfedge()->twin();
    HalfedgeRef he3_fast = e->halfedge()->twin();
    do {
        he3 = he3->next();
        he3_fast = he3_fast->next()->next();
    } while (he3_fast->id() != e->halfedge()->twin()->id() && he3_fast->next()->id() != e->halfedge()->twin()->id());
    HalfedgeRef he4 = he3->next();
    std::cout << "he4 id: " + std::to_string(he4->id()) + "\n";

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
    // he2->vertex()->_halfedge = new_he;

    he3->next() = new_he_twin;
    new_he_twin->next() = he2;
    new_he_twin->vertex() = he4->vertex();
    new_he_twin->vertex()->_halfedge = new_he_twin;
    // he3->vertex()->_halfedge = new_he_twin;

    he5->next() = e->halfedge()->twin()->next();
    he6->next() = e->halfedge()->next();
    e->halfedge()->vertex()->_halfedge = he5->next();
    e->halfedge()->twin()->vertex()->_halfedge = he6->next();
    
    std::cout << "HERE c\n";
    
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

    std::cout << "New edge + face info:\n";
    curr = new_e->halfedge();
    do {
        std::cout << "he: " + std::to_string(curr->id()) + " v: " + std::to_string(curr->vertex()->id()) + " (" + std::to_string(curr->vertex()->halfedge()->id()) + " -> " + std::to_string(curr->vertex()->halfedge()->vertex()->id()) + ") f: " + std::to_string(curr->face()->id()) + "\n";
        curr = curr->next();
    } while (curr->id() != new_e->halfedge()->id());

    std::cout << "New twin edge + face info:\n";
    curr = new_e->halfedge()->twin();
    std::cout << "he: " + std::to_string(curr->id()) + "\n";
    std::cout << "v: " + std::to_string(curr->vertex()->id()) + "\n";
    std::cout << "f: " + std::to_string(curr->face()->id()) + "\n";
    do {
        std::cout << "he: " + std::to_string(curr->id()) + " v: " + std::to_string(curr->vertex()->id()) + " (" + std::to_string(curr->vertex()->halfedge()->id()) + " -> " + std::to_string(curr->vertex()->halfedge()->vertex()->id()) + ") f: " + std::to_string(curr->face()->id()) + "\n";
        curr = curr->next();
    } while (curr->id() != new_e->halfedge()->twin()->id());

    std::cout << std::to_string(e->halfedge()->id()) + " points to " + std::to_string(e->halfedge()->next()->id()) + "\n";
    std::cout << std::to_string(e->halfedge()->twin()->id()) + " points to " + std::to_string(e->halfedge()->twin()->next()->id()) + "\n";

    // Remove the original edge.
    std::cout << "eraseing edge: " + std::to_string(e->id()) + "\n";
    std::cout << "eraseing he: " + std::to_string(e->halfedge()->id()) + "\n";
    std::cout << "eraseing he: " + std::to_string(e->halfedge()->twin()->id()) + "\n";
    erase(e);
    erase(e->halfedge());
    erase(e->halfedge()->twin());
    validate(); // DEBUG
    std::cout << "HERE e\n";
    // Return a pointer to the new flipped edge.
    return new_e;
    // return std::nullopt;
}

/*
    This method should split the given edge and return an iterator to the
    newly inserted vertex. The halfedge of this vertex should point along
    the edge that was split, rather than the new edges.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::split_edge(Halfedge_Mesh::EdgeRef e) {

    // (void)e;
    // NOTE: Getting face/vertex halfedge don't point to face/vertex errors, but 
    // print statements suggest otherwise.
    // Create three new edges and their halfedges, a new vertex, and two new faces.
    // We pretend that e->halfedge points "north".
    EdgeRef west_e = new_edge();
    EdgeRef south_e = new_edge();
    EdgeRef east_e = new_edge();
    HalfedgeRef west_he = new_halfedge();
    HalfedgeRef west_he_twin = new_halfedge();
    HalfedgeRef south_he = new_halfedge();
    HalfedgeRef south_he_twin = new_halfedge();
    HalfedgeRef east_he = new_halfedge();
    HalfedgeRef east_he_twin = new_halfedge();
    VertexRef v = new_vertex();
    FaceRef bl_f = new_face();
    FaceRef br_f = new_face();
    VertexRef north_v = e->halfedge()->twin()->vertex();
    VertexRef west_v = e->halfedge()->next()->next()->vertex();
    VertexRef south_v = e->halfedge()->vertex();
    VertexRef east_v = e->halfedge()->twin()->next()->next()->vertex();
    HalfedgeRef ul_he = e->halfedge()->next();
    HalfedgeRef ll_he = e->halfedge()->next()->next();
    HalfedgeRef lr_he = e->halfedge()->twin()->next();
    HalfedgeRef ur_he = e->halfedge()->twin()->next()->next();

    // Reassign the west, south, and east edges and halfedges.
    west_he->_next = ll_he;
    west_he->_twin = west_he_twin;
    west_he->_vertex = v;
    west_he->_edge = west_e;
    west_he->_face = bl_f;
    west_he_twin->_next = e->halfedge();
    west_he_twin->_twin = west_he;
    west_he_twin->_vertex = west_v;
    west_he_twin->_edge = west_e;
    west_he_twin->_face = e->halfedge()->face();
    west_e->_halfedge = west_he;

    south_he->_next = lr_he;
    south_he->_twin = south_he_twin;
    south_he->_vertex = v;
    south_he->_edge = south_e;
    south_he->_face = br_f;
    south_he_twin->_next = west_he;
    south_he_twin->_twin = south_he;
    south_he_twin->_vertex = south_v;
    south_he_twin->_edge = south_e;
    south_he_twin->_face = bl_f;
    south_e->_halfedge = south_he;

    east_he->_next = ur_he;
    east_he->_twin = east_he_twin;
    east_he->_vertex = v;
    east_he->_edge = east_e;
    east_he->_face = e->halfedge()->twin()->face();
    east_he_twin->_next = south_he;
    east_he_twin->_twin = east_he;
    east_he_twin->_vertex = east_v;
    east_he_twin->_edge = east_e;
    east_he_twin->_face = br_f;
    east_e->_halfedge = east_he;

    // Reassign the bottom-left and bottom-right faces. We can keep the same halfedge
    // for the top-left and top-right faces.
    bl_f->_halfedge = west_he;
    br_f->_halfedge = south_he;
    e->halfedge()->twin()->face()->_halfedge = east_he;

    // Reassign the new vertex.
    v->_halfedge = e->halfedge();

    // Update boundary halfedge next pointers.
    ll_he->_next = south_he_twin;
    ul_he->_next = west_he_twin;
    lr_he->_next = east_he_twin;

    // We will use the original edge and halfedges that we had for the north edge and halfedges.
    // We only need to udpate the twin halfedge because e->halfedge points towards
    // unchanged regions.
    e->halfedge()->_vertex = v;
    e->halfedge()->twin()->_next = east_he;

    // HALF EDGE NEXT CHECKER.
    HalfedgeRef curr = e->halfedge();
    std::cout << "\nnorth he " + std::to_string(curr->id()) + " -> he " + std::to_string(curr->next()->id()) + "\n";
    curr = curr->next();
    std::cout << "northwest he " + std::to_string(curr->id()) + " -> he " + std::to_string(curr->next()->id()) + "\n";
    curr = curr->next();
    std::cout << "west he twin " + std::to_string(curr->id()) + " -> he " + std::to_string(curr->next()->id()) + "\n";

    curr = west_he;
    std::cout << "\nwest he " + std::to_string(curr->id()) + " -> he " + std::to_string(curr->next()->id()) + "\n";
    curr = curr->next();
    std::cout << "southwest he " + std::to_string(curr->id()) + " -> he " + std::to_string(curr->next()->id()) + "\n";
    curr = curr->next();
    std::cout << "south he twin " + std::to_string(curr->id()) + " -> he " + std::to_string(curr->next()->id()) + "\n";

    curr = south_he;
    std::cout << "\nsouth he " + std::to_string(curr->id()) + " -> he " + std::to_string(curr->next()->id()) + "\n";
    curr = curr->next();
    std::cout << "southeast he " + std::to_string(curr->id()) + " -> he " + std::to_string(curr->next()->id()) + "\n";
    curr = curr->next();
    std::cout << "east he twin " + std::to_string(curr->id()) + " -> he " + std::to_string(curr->next()->id()) + "\n";

    curr = east_he;
    std::cout << "\neast he " + std::to_string(curr->id()) + " -> he " + std::to_string(curr->next()->id()) + "\n";
    curr = curr->next();
    std::cout << "northeast he " + std::to_string(curr->id()) + " -> he " + std::to_string(curr->next()->id()) + "\n";
    curr = curr->next();
    std::cout << "north he twin " + std::to_string(curr->id()) + " -> he " + std::to_string(curr->next()->id()) + "\n";

    // FACE CHECKER.
    FaceRef faceref = e->halfedge()->face();
    std::cout << "\ntl face " + std::to_string(faceref->id()) + " -> " + std::to_string(faceref->halfedge()->face()->id()) + "\n";
    faceref = east_he->face();
    std::cout << "tr face " + std::to_string(faceref->id()) + " -> " + std::to_string(faceref->halfedge()->face()->id()) + "\n";
    faceref = south_he->face();
    std::cout << "lr face " + std::to_string(faceref->id()) + " -> " + std::to_string(faceref->halfedge()->face()->id()) + "\n";
    faceref = west_he->face();
    std::cout << "ll face " + std::to_string(faceref->id()) + " -> " + std::to_string(faceref->halfedge()->face()->id()) + "\n";

    // VERTEX CHECKER.
    VertexRef vref = v;
    std::cout << "\ncenter v " + std::to_string(vref->id()) + " -> " + std::to_string(vref->halfedge()->vertex()->id()) + "\n";
    vref = north_v;
    std::cout << "north v " + std::to_string(vref->id()) + " -> " + std::to_string(vref->halfedge()->vertex()->id()) + "\n";
    vref = west_v;
    std::cout << "west v " + std::to_string(vref->id()) + " -> " + std::to_string(vref->halfedge()->vertex()->id()) + "\n";
    vref = south_v;
    std::cout << "south v " + std::to_string(vref->id()) + " -> " + std::to_string(vref->halfedge()->vertex()->id()) + "\n";
    vref = east_v;
    std::cout << "east v " + std::to_string(vref->id()) + " -> " + std::to_string(vref->halfedge()->vertex()->id()) + "\n";

    // Return a pointer to the vertex.
    validate();
    return v;
    // return std::nullopt;
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

    (void)f;
    return std::nullopt;
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
}

/*
    Splits all non-triangular faces into triangles.
*/
void Halfedge_Mesh::triangulate() {
    // TODO: Check error "A face has a halfedge which does not point to that face!"ome
    // Iterate and store pointers to all original mesh faces.
    std::vector<FaceRef> og_faces;
    FaceRef curr_face = faces_begin();
    do {
        std::cout << "face: " + std::to_string(curr_face->id()) + "\n";
        og_faces.push_back(curr_face);
        curr_face++;
    } while (curr_face != faces_end());

    // For each face, choose a vertex u. Iterate through all other vertices v around
    // the face and add an edge u -> v if there doesn't already exist one.
    for (int f = 0; f < (int)og_faces.size(); f++) {
        curr_face = og_faces[f];
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
        face_hes[face_hes.size() - 1]->_next = new_halfedges[new_halfedges.size() - 1];

        // Reassign the halfedge of the original face.
        curr_face->_halfedge = new_halfedges[new_halfedges.size() - 1];
    }
    validate();

    // CHECK FACE HALFEDGE POINTERS.
    for (FaceRef f = faces_begin(); f != faces_end(); f++) {
        std::cout << "f " + std::to_string(f->id()) + " -> he " + std::to_string(f->halfedge()->id()) + " -> f " + std::to_string(f->halfedge()->face()->id()) + "\n";
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
        VertexRef u = e->halfedge()->vertex();
        VertexRef v = e->halfedge()->next()->vertex();
        Vec3 midpoint = (u->pos + v->pos) / 2;
        e->new_pos = midpoint;
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

    // Edges

    // Vertices
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
        // -> Build the 3x3 linear system whose solution minimizes the quadric error
        //    associated with these two endpoints.
        // -> Use this system to solve for the optimal position, and store it in
        //    Edge_Record::optimal.
        // -> Also store the cost associated with collapsing this edge in
        //    Edge_Record::cost.
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
    // -> Compute an initial quadric for each vertex as the sum of the quadrics
    //    associated with the incident faces, storing it in vertex_quadrics
    // -> Build a priority queue of edges according to their quadric error cost,
    //    i.e., by building an Edge_Record for each edge and sticking it in the
    //    queue. You may want to use the above PQueue<Edge_Record> for this.
    // -> Until we reach the target edge budget, collapse the best edge. Remember
    //    to remove from the queue any edge that touches the collapsing edge
    //    BEFORE it gets collapsed, and add back into the queue any edge touching
    //    the collapsed vertex AFTER it's been collapsed. Also remember to assign
    //    a quadric to the collapsed vertex, and to pop the collapsed edge off the
    //    top of the queue.

    // Note: if you erase elements in a local operation, they will not be actually deleted
    // until do_erase or validate are called. This is to facilitate checking
    // for dangling references to elements that will be erased.
    // The rest of the codebase will automatically call validate() after each op,
    // but here simply calling collapse_edge() will not erase the elements.
    // You should use collapse_edge_erase() instead for the desired behavior.

    return false;
}
