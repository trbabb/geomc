#include <limits>
#include <vector>

#include <geomc/linalg/Vec.h>


namespace geom {

namespace detail {


// returns 0 if unequal, 1 if same orientation, -1 if opposite orientation
// pass a mutable simplex to s1 plz, kthx
int simplex_comparison(const index_t* s0, index_t* s1, index_t n) {
    // reorder s1 to look like s0, keeping track of how many transpositions are made
    // to determine parity.
    int swaps = 0;
    for (index_t i = 0; i < n; i++) {
        // find s0[i] in s1 and swap positions.
        for (index_t j = i; j < n; j++) {
            if (s0[i] == s1[j]) {
                std::swap(s1[i], s1[j]);
                swaps++;
                goto FOUNDIT;
            }
        }
        return 0; // got to end of list without finding current element.
        FOUNDIT:
    }
    return (swaps & 1) ? -1 : 1;
}


template <typename T, index_t N>
struct Edge {
    index_t v[N-1]; // in ccw order. 
                    // we use indexes instead of pointers because resizing a std::vector
                    // invalidates pointers and iterators.
};


template <typename T, index_t N>
struct Face {
    index_t v[N]; // in ccw order.
    Vec<T,N> n;   // face normal.
};


} // namespace detail


// todo: simplex_comparison can be a hard equality check for n=2
// todo: for n=2, we know which edges (verts) are to be patched exactly.
//       there are two of them, and dynamic list is unnecessary.
// todo: does this algorithm REALLY generalize this way?
// todo: xxx: patch up gjk_intersect() to return its internal simplex
// todo: write an optimized N=2 implementation and then verify it against the general one

// observation: The number of faces is strictly increasing.
//              you don't need to delete a face (O(n)); they can be clobbered.
//              use a skip list.
// note that there is no reason to extend this algorithm to find a
// separation axis, because GJK already gives that to us.


template <typename T, index_t N>
bool separation_axis(const Convex<T,N>& shape_a,
                     const Convex<T,N>& shape_b,
                     Vec<T,N>* overlap_axis,
                     double fractional_tolerance = 0.001,
                     index_t iteration_limit = -1) {
    Simplex<T,N> splex;
    
    // assignment below is intended.
    if (not gjk_intersect(shape_a, shape_b, overlap_axis, &splex)) {
        return false;
    }
    
    std::vector< Vec<T,N> >          verts(splex.pts, splex.pts + splex.n);
    std::vector< detail::Edge<T,N> > edges;
    std::vector< detail::Face<T,N> > faces;
    
    Vec<T,N> splex_vtx_buf[N-1];
    
    // initialize the `faces` vec with those of the GJK simplex.
    for (index_t i = 0; i <= N; i++) {
        // choose a vertex to exclude, making a face, and construct it.
        Vec<T,N> A = verts[i];
        Face<T,N> face;
        for (index_t j = 0; j < N; j++) {
            face.v[j] = (i + j + 1) % N;
            if (j > 0) splex_vtx_buf[j-1] = verts[face.v[0]] - verts[i];
        }
        
        // correct the winding.
        Vec<T,N> n = orthogonal(splex_vtx_buf);
        if (n.dot(A - verts[face.v[0]]) > 0) {
            // winding is backwards; normal points inwards. 
            // make an arbitrary transposition to flip it.
            std::swap(face.v[0], face.v[1]);
            n = -n;
        }
        face.n = n.unit();
        faces.push_back(face);
    }
    
    std::vector< detail::Face<T,N> >::iterator best_face;
    Vec<T,N> last_proj;
    bool looped = false;
    
    // iterate on the face list, choosing the closest one and expanding it
    while (index_t k = 0; k < iteration_limit or iteration_limit < 1; k++) {
        
        // find the face closest to the origin.
        T d = std::numeric_limits<T>::max();
        for (auto f = faces.begin(); f != faces.end(); ++f) {
            Vec<T,N> v = verts[f->v[0]]; // a vertex on the face
            T d_this = std::abs(f->n.dot(v));
            if (d_this < d) {
                d = d_this;
                best_face = f;
            }
        }
        
        // project the origin to the closest face.
        *overlap_axis = verts[best_face->v[0]].projectOn(best_face->n);
        
        // if our estimated closest point is not sufficiently different from
        // our last estimate, decide that we've converged and quit.
        if ((looped and fractional_tolerance > 0 and 
                (last_proj - overlap_axis).mag() / overlap_axis.mag() < fractional_tolerance)
                or *overlap_axis == Vec<T,N>::zeros) {
            break;
        }
        
        last_proj = *overlap_axis;
        
        // look in a direction normal to that face
        // and find a new point on the minkowski difference.
        Vec<T,N> minkowski_pt = shape_a.convex_support( best_face->n) - 
                                shape_b.convex_support(-best_face->n);
        
        // finish if we add a duplicate point (possible for polytopes)
        // such a point would produce a degenerate face, i.e. we aren't going
        // to find any more new faces. this is a linear algorithm, but I am 
        // not sure whether a hashset construction/test is actually worth it.
        for (auto v : verts) if (v == minkowski_pt) break;
        
        verts.push_back(minkowski_pt);
        
        // now delete all the faces that fall "behind" the new point.
        for (auto f = faces.begin(); f != faces.end();) {
            // if this face is below the new point...
            if (f->n.dot(minkowski_pt - verts[f->pts[0]]) > 0) {
                // add all the edges of this face to the boundary edge list.
                // (all the non-internal edges are hole boundary, which must each
                //  become part of a new face).
                for (index_t i = 0; i < N; i++) {
                    // find a vertex to exclude from the face, making an edge.
                    Edge<T,N> e;
                    for (index_t j = 0; j < N - 1; j++) {
                        e.v[j] = f->pts[(i + j + 1) % N];
                        if (N > 3 and j > 0) {
                            splex_vtx_buf[j-1] = verts[e.v[0]] - verts[e.v[j]];
                        }
                    }
                    
                    // correct the edge winding, if necessary.
                    if (N > 3) {
                        splex_vtx_buf[N-2] = f->n;
                        // the parity here is arbitrary; all that matters is that it's consistent.
                        Vec<T,N> n_edge = orthogonal(splex_vtx_buf); 
                        if (n_edge.dot(verts[e.v[0]] - splex_vtx_buf[0]) > 0) {
                            std::swap(e.v[0], e.v[1]);
                        }
                    }
                    
                    // add the edge.
                    // if there exists an edge in the list with opposite winding, they "annihilate"
                    // we should actually never see a collision with equal winding, so if we do, that's bad.
                    Edge<T,N> e_tmp;
                    bool annihilated = false;
                    for (auto e_other = edges.begin(); e_other != edges.end(); ++e_other) {
                        e_tmp = *e_other;
                        if (simplex_comparison(e, e_tmp, N-1) == -1) {
                            edges.erase(e_other);
                            annihilated = true;
                            break;
                        }
                    }
                    if (not annihilated) edges.push_back(e);
                }
                // delete the occluded face.
                f = faces.erase(f);
            } else {
                ++f;
            }
        }
        
        // patch up the hole by constructing new faces,
        // connecting each hole-adjacent edge to the new vertex.
        index_t new_vert = verts.size() - 1;
        for (auto e = edges.begin(); e != edges.end(); ++e) {
            Face<T,N> f;
            std::copy(e->v, e->v + N - 1, f.v);
            f.v[N-1] = new_vert;
            
            // compute face normal
            for (index_t i = 0; i < N - 1; i++) {
                splex_vtx_buf[i] = verts[face[i]] - minkowski_pt;
            }
            f.n = orthogonal(splex_vtx_buf).unit();
            
            // correct the face winding, if necessary.
            if (N > 3) {
                // we know the origin is inside the polytope, so we will
                // use that to test whether the normal is pointing the right way.
                // so here we treat minkowski_pt as a vector pointing from inside
                // the polytope to a point on the face simplex.
                if (f.n.dot(minkowski_pt) < 0) {
                    f.n = -f.n;
                    std::swap(f.v[0], f.v[1]);
                }
            }
            
            faces.push_back(f);
        }
        
        edges.clear();
        looped = true;
    }
    
    return true;
}
    
    
    
}