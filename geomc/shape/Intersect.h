/* 
 * File:   Intersect.h
 * Author: tbabb
 *
 * Created on October 14, 2014, 11:30 PM
 */

#ifndef INTERSECT_H
#define	INTERSECT_H

#include <limits.h>
#include <geomc/linalg/Vec.h>
#include <geomc/linalg/Orthogonal.h>
#include <geomc/shape/Bounded.h>

#include <vector>
#include <map>
 
// xxx debug
#define EMIT(...) if(emit_debug) emit(__VA_ARGS__)

// todo: can we get faster/more stable/simpler results by solving for barycentric coords?

namespace geom {
    
namespace detail {
    
    
    template <typename T, index_t N>
    struct Simplex {
        Simplex():n(0) {}
        
        Vec<T,N> pts[N+1];
        index_t  n;
        
        inline void insert(const Vec<T,N> &p) {
            if (n <= N) pts[n++] = p;
        }
        
        inline void remove(index_t i) {
            for (index_t c = i; c < n - 1; c++) {
                pts[c] = pts[c+1];
            }
            n -= 1;
        }
        
        inline const Vec<T,N>& operator[](index_t i) const {
            return pts[i];
        }
        
        inline Vec<T,N>& operator[](index_t i) {
            return pts[i];
        }
    };
    
    
    // Store a sub-simplex, keeping track of its null basis. 
    // The sub-simplex implicitly includes an unstored vertex A, which in GJK
    // is common to all tested sub-simplexes. 
    template <typename T, index_t N>
    struct SubSimplex {
        
        SubSimplex():n(0) {}
        
        const Vec<T,N> *pts[N+1]; // array of <= N pointers to const vec
        Vec<T,N> null_basis[N];   // basis for orthogonal complement of spanned space
        index_t n; // number of stored points (i.e., excludes "A"). equal to the dimension spanned by `this`.
        
        // create a new, complete simplex from `s`. The final vertex is of `s`
        // is "A", and is not stored. 
        SubSimplex& createFrom(const Simplex<T,N> &s) {
            // the last vertex is assumed to be `A`
            n = s.n - 1;
            for (index_t i = 0; i < n; i++) {
                pts[i] = &s[i];
            }
            // compute null space (if it's not one of {R^N, 0})
            if (n > 0 and n < N) {
                // compute the null space of this simplex.
                Vec<T,N> A = s[s.n - 1];
                Vec<T,N> *simplex_basis = null_basis + N - n;
                for (index_t i = 0; i < n; i++) {
                    simplex_basis[i] = s[i] - A;
                }
                nullspace(simplex_basis, n, null_basis);
            }
            return *this;
        }
        
        // create a new sub-simplex of lower dimension by excluding the vertex with
        // index `excluded`. The point `A` is taken to be the final vertex of this
        // sub simplex, and is used in calculating the null space and normal.
        SubSimplex& createFrom(const SubSimplex &s, index_t excluded, const Vec<T,N> &A) {
            n = s.n - 1;
            index_t j = 0;
            for (index_t i = 0; i < s.n; i++) {
                if (i == excluded) continue;
                pts[j++] = s.pts[i];
            }
            // compute the simplex normal
            std::copy(s.null_basis, s.null_basis + N - s.n, null_basis);
            Vec<T,N> *simplex_basis = null_basis + N - s.n; // high part of null_basis array
            for (index_t i = 0; i < n; i++) {
                simplex_basis[i] = *(pts[i]) - A;
            }
            if (n > 0) {
                // the normal lies orthogonal to the null space of the parent
                // and orthogonal to the space of the child. null_space now
                // contains both.
                Vec<T,N> normal = orthogonal(null_basis);
                // flip the normal if it points "inside"
                if (normal.dot(s[excluded] - A) > 0) normal = -normal;
                // extend the null basis.
                simplex_basis[0] = normal;
            } else {
                // we are a point; parent simplex is a line; there is only one normal
                simplex_basis[0] = A - s[0];
            }
            return *this;
        }
        
        const Vec<T,N>& operator[](index_t i) const {
            return *(pts[i]);
        }
        
        inline Vec<T,N> getNormal() const {
            return null_basis[N-n-1];
        }
        
    };
    
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
            FOUNDIT: continue;
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
        
        inline bool operator==(const Face<T,N>& other) {
            for (index_t i = 0; i < N; i++) if (v[i] != other.v[i]) return false;
            return true;
        }
    };
    
    
    // statically-sized/allocated queue
    // (behavior undefined if size grows > N)
    template <typename T, index_t N>
    class StaticQueue {
        T q[N];
        index_t front; // points at first item
        index_t back;  // points at empty slot
        
        public:
        
        index_t size;
        
        StaticQueue():front(0),back(0),size(0) {}
        
        inline T* create_back() {
            T *ret = q + back;
            back = (back + 1) % N;
            size++;
            return ret;
        }
        
        inline T& peek_front() {
            return q[front];
        }
        
        inline void destroy_front() {
            front = (front + 1) % N;
            size--;
        }
        
        inline void destroy_back() {
            back = positive_mod(back - 1, N);
            size--;
        }
    };
    
    
    constexpr index_t gjk_queue_size_fac(index_t n) {
        return (n < 2) ? 1 : n * gjk_queue_size_fac(n-1);
    }
    
    
    template <typename T, index_t N>
    bool gjk_direction_to_origin(detail::SubSimplex<T,N> &boundary, const Vec<T,N> &A, Vec<T,N> *d) {
        // find the direction to the origin. 
        // assumes a simplex which has already passed containment tests.
        // project the direction to the origin onto the null space of our simplex.
        // the final `else` would be sufficient by itself, but in some cases
        // we can be clever and save some redundant calculations.
        index_t null_dim = N - boundary.n;
        if (null_dim == N) {
            // single point case; "projection" is trivial.
            *d = -A;
        } else if (null_dim == 1) {
            // the simplex spans a hyperplane; we already have the normal.
            *d = boundary.getNormal();
            // if we created the plane by removing a vertex, we already know the 
            // normal is facing toward the origin. otherwise, we need to check 
            // and flip it if necessary.
            if (d->dot(A) > 0) *d *= -1;
        } else if (null_dim == 0) {
            // simplex forms a complete basis; no vertexes were eliminated.
            // in other words, no plane test failed. The origin is inside
            // this simplex!
            return true;
        } else {
            // explicitly project the direction to the origin onto the null space
            // of the boundary simplex.
            if (N != 3) {
                *d = Vec<T,N>::zeros;
                
                // either project onto the null basis directly, or subtract a vector
                // orthogonal to the null basis from -A, according to whichever 
                // projection is simpler.
                bool swap_proj = false;
                index_t proj_n = null_dim;
                Vec<T,N> *proj_basis = boundary.null_basis;

                if ( null_dim > boundary.n ) {
                    proj_basis = boundary.null_basis + N - boundary.n;
                    proj_n = boundary.n;
                    swap_proj = true;
                    for (index_t i = 0; i < boundary.n; i++) {
                        proj_basis[i] = boundary[i] - A;
                    }
                }

                // gram-schmidt orthonormalization.
                for (index_t b = 0; b < proj_n; b++) {
                    for (index_t i = 0; i < b; i++) {
                        proj_basis[b] -= proj_basis[b].projectOn(proj_basis[i]);
                    }
                    *d += -A.projectOn(proj_basis[b]);
                }

                if (swap_proj) *d = -A - *d;
            } else {
                // cross product is faster/more numerically stable in 3D, but
                // results are in principle identical. 
                
                // Vec<T,N> b = A - boundary[0];
                // Vec<T,N> c = b ^ -A;
                // *d = c ^ b; 
                
                // We use orthogonal() instead of ^ because this is a 
                // dimension-agnostic function. orthogonal inlines to ^ anyway.
                Vec<T,N> vs[3] = {(T)0, A - boundary[0], -A};
                vs[0] = orthogonal(vs+1);
                *d = orthogonal(vs);
            }
        }
        return false;
    }
    
    
    /* todo:
     * The only vertex that might be a boundary vertex is A. before we know A is the 
     * boundary vertex, we must check that the origin does not belong to an edge.
     * There are at most only N edges in the entire simplex. However, there are potentially
     * an absurd number of triangles, and we will check as many as two times that number of edges.
     * This is obviously very redundant because the null space of an edge (and uniquely so
     * for edges only) is not dependent on the parent, since there is only one degree 
     * of freedom. However, it is not clear to me how to avoid redundant checks in a way
     * that is not O(num children of triangles) anyway. Since a plane test is so damn
     * cheap, it may not be a big deal. (Also note that each edge needs a complete null basis,
     * which would have to be constructed for each edge). At best, this might save
     * us some storage.
     */
    
    // todo: check for degenerate simplexes and ignore them.
    
    template <typename T, index_t N>
    bool gjk_simplex_nearest_origin(detail::Simplex<T,N> *simplex, Vec<T,N> *d) {
        typedef detail::SubSimplex<T,N> SS;
        // queue for sub-simplexes to check.
        // (sub-simplexes store their null basis, and ptrs to their verts).
        detail::StaticQueue<SS, detail::gjk_queue_size_fac(N) + 1> q;

        Vec<T,N> A = simplex->pts[simplex->n - 1];
        SS boundary;
        
        // enqueue this simplex, and its null basis 
        q.create_back()->createFrom(*simplex);

        // process the "frontier"
        while (q.size > 0) {
            SS &s = q.peek_front(); // note that `s` includes `A` implicitly.
            // a simplex is a "boundary" simplex if the origin projects down its normal within its bounds.
            // if any of the simplex's edge planes do not contain the origin, this cannot be true.
            // however, that edge itself may be a boundary simplex.
            bool parent_is_boundary = true;
            for (index_t i = 0; i < s.n; i++) {
                SS &child = q.create_back()->createFrom(s, i, A);
                Vec<T,N> normal = child.getNormal();

                if (A.dot(normal) < 0) { 
                    // the origin is "outside" this sub-simplex, and therefore it
                    // may be (or contain) a boundary simplex. We also know the parent
                    // cannot be the boundary simplex. 
                    parent_is_boundary = false;
                } else {
                    // no need to check the grandchildren.
                    q.destroy_back();
                }
            }
            if (parent_is_boundary) {
                boundary = s;
                break;
            } else {
                // this simplex is checked; pop it off the stack.
                q.destroy_front();
            }
        }
        
        bool containment = gjk_direction_to_origin(boundary, A, d);
        
        // put the new simplex back into the return variable.
        // we use a temp because the boundary actually points into `simplex`.
        Simplex<T,N> tmp;
        tmp.n = boundary.n + 1;
        tmp.pts[boundary.n] = A;
        for (index_t i = 0; i < boundary.n; i++) {
            tmp[i] = boundary[i];
        }
        *simplex = tmp;
        
        return containment;
    }
    
    template <typename T, index_t N>
    struct UnstructuredPointcloud : virtual public Convex<T,N> {
        const Vec<T,N> *pts;
        index_t n;
        
        UnstructuredPointcloud(const Vec<T,N> *pts, index_t n):pts(pts),n(n) {}
        
        Vec<T,N> convexSupport(Vec<T,N> d) const {
            T largest_dot = std::numeric_limits<T>::lowest();
            Vec<T,N> pt;
            for (index_t i = 0; i < n; i++) {
                T a = d.dot(pts[i]);
                if (a > largest_dot) {
                    largest_dot = a;
                    pt = pts[i];
                }
            }
            return pt;
        }
    };
    

    template <typename T, index_t N>
    bool gjk_intersect(const Convex<T,N> &shape_a, 
                       const Convex<T,N> &shape_b,
                       Vec<T,N> *overlap_axis,
                       detail::Simplex<T,N> *s) {
        // choose an arbitrary initial direction.
        Vec<T,N> initial;
        if (overlap_axis and *overlap_axis != Vec<T,N>::zeros) 
            initial = *overlap_axis;
        else 
            initial[0] = 1;
        
        Vec<T,N> a = shape_a.convexSupport( initial) - 
                     shape_b.convexSupport(-initial);
        Vec<T,N> d = -a;
        s->n = 0;
        s->insert(a);
        
        // supported by empirical measurement, plus a safe margin:
        const index_t cutoff = 10 * (1 << (N - 2)); 
        index_t i = 0;
        
        while (true) {
            a = shape_a.convexSupport( d) - 
                shape_b.convexSupport(-d);
            if (a.dot(d) < 0) {
                if (overlap_axis) *overlap_axis = -a.projectOn(d);
                return false;
            }
            s->insert(a);
            
            if (detail::gjk_simplex_nearest_origin(s, &d)) {
                if (overlap_axis) *overlap_axis = d;
                return true;
            }
            if (++i > cutoff) { return true; }
        }
    }
    

} // namespace detail
 

// todo: what if I want to do stuff like dilate or erode the minkowski difference?
//       what if I want to do stuff to only one of the shapes?

/**
 * @ingroup shape
 * Use the Gilbert-Johnson-Keerthi algorithm to determine if two convex shapes
 * overlap, and return an axis of overlap.
 * @param shape_a A convex shape.
 * @param shape_b A convex shape.
 * @param overlap_axis An initial axis along which to test for overlap; upon completion
 * this will be populated with an axis (not necessarily minimal) of overlap or separation.
 * Pass null to choose an arbitrary initial axis.
 * @return `true` if and only if `shape_a` and `shape_b` overlap.
 */
template <typename T, index_t N>
#ifndef PARSING_DOXYGEN
inline 
#endif 
bool gjk_intersect(const Convex<T,N> &shape_a, 
                   const Convex<T,N> &shape_b,
                   Vec<T,N> *overlap_axis=NULL) {
    detail::Simplex<T,N> s;
    return detail::gjk_intersect(shape_a, shape_b, overlap_axis, &s);
}


template <typename T, index_t N>
inline bool gjk_intersect(const Vec<T,N> *pts_a, index_t n_a,
                          const Vec<T,N> *pts_b, index_t n_b,
                          Vec<T,N> *overlap_axis=NULL) {
    detail::UnstructuredPointcloud<T,N> a(pts_a, n_a);
    detail::UnstructuredPointcloud<T,N> b(pts_b, n_b);
    return gjk_intersect(a, b, overlap_axis);
}

// todo: write an optimized N=2 implementation and then verify it against the general one
//       - edge comparison is a hard equality check
//       - there can only ever be exactly two unpatched hole edges.
// todo: optimize face insertion
// todo: allow buffer containers to be passed.
// todo: test at what point hashtable construction is faster for duplicate pt tests.
// todo: keep the faces in sorted order by distance. insertion is O(lg(n)) and
//       searching is O(1). std::list means no shuffling.

// observation: The number of faces is strictly increasing.
//              you don't need to delete a face [O(n)]; they can be clobbered.
//              use a ~~skip list~~ std::list.
// note that there is no reason to extend this algorithm to find a
// separation axis, because GJK already gives that to us.

// begin debug {

template <typename T, index_t N>
void emit(const detail::Face<T,N>& f, bool highlight=false) {
    if (highlight) std::cout << "h ";
    std::cout << "f ";
    for (index_t i = 0; i < N; i++) {
        std::cout << f.v[i] << " ";
    }
    std::cout << "n " << f.n << "\n";
}

template <typename T, index_t N>
void emit(const detail::Edge<T,N>& e, bool highlight=false) {
    if (highlight) std::cout << "h ";
    std::cout << "e ";
    for (index_t i = 0; i < N-1; i++) {
        std::cout << e.v[i] << " ";
    }
    std::cout << "\n";
}

template <typename T, index_t N>
void emit(const Vec<T,N>& v) {
    std::cout << "v ";
    for (index_t i = 0; i < N; i++) std::cout << v[i] << " ";
    std::cout << "\n";
}

// } end debug


bool disjoint_separation_axis(const Convex<T,N>& shape_a,
                              const Convex<T,N>& shape_b, 
                              Vec<T,N> *overlap_axis,
                              const detail::Simplex<T,N>& splex,
                              double fractional_tolerance = 0.001,
                              index_t iteration_limit = -1) {
    std::vector< Vec<T,N> >                  verts;
    std::vector< detail::Edge<T,N> >         edges; // problem. subsimplexes have ptrs to verts, but a std::vector's ptr may be invalidated. :(
    std::multimap< T, detail::SubSimplex<T,N> > faces;
}


/**
 * @ingroup shape
 * Use the Expanding Polytope Algorithm to find a minimum translation vector that would bring shape B into contact with A.
 * 
 * @param shape_a A convex shape.
 * @param shape_b A convex shape.
 * @param overlap_axis An initial test axis, and return variable for the separation vector.
 * @param fractional_tolerance Convergence condition: When the magnitude of the difference
 * between iterations is less than this fraction of the total separation, terminate the algorithm.
 * @param iteration_limit Hard limit on number of convergence iterations. Pass -1 (default) to 
 * impose no hard limit.
 * @return `true` if and only if the two shapes overlap.
 */
template <typename T, index_t N>
bool minimal_separation_axis(const Convex<T,N>& shape_a,
                             const Convex<T,N>& shape_b,
                             Vec<T,N>* overlap_axis,
                             double fractional_tolerance = 0.001,
                             index_t iteration_limit = -1,
                             bool emit_debug=false) { // xxx debug
    detail::Simplex<T,N> splex;
    
    if (not detail::gjk_intersect(shape_a, shape_b, overlap_axis, &splex)) {
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
        
            EMIT(A);
        
        detail::Face<T,N> face;
        for (index_t j = 0; j < N; j++) {
            face.v[j] = (i + j + 1) % (N + 1);
            if (j > 0) splex_vtx_buf[j-1] = verts[face.v[0]] - verts[face.v[j]];
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
        
            EMIT(face);
    }
    Vec<T,N> last_proj;
    bool looped = false;
    
    // iterate on the face list, choosing the closest one and expanding it
    for (index_t k = 0; iteration_limit < 1 or k < iteration_limit; k++) {
        typename std::vector< detail::Face<T,N> >::iterator best_face;
        
        if (emit_debug) std::cout << "=\n";
        
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
        
            // debug: emit current polytope, highlighting the closest face
            for (auto f = faces.begin(); f != faces.end(); ++f) EMIT(*f, f == best_face);
        
        // project the origin to the closest face.
        *overlap_axis = verts[best_face->v[0]].projectOn(best_face->n);
        
        // if our estimated closest point is not sufficiently different from
        // our last estimate (or if the origin is on the current face),
        // decide that we've converged and quit.
        if ((looped and fractional_tolerance > 0 and 
                (last_proj - *overlap_axis).mag() / overlap_axis->mag() < fractional_tolerance) 
                // origin is exactly on the face:
                or *overlap_axis == Vec<T,N>::zeros) {
            break;
        }
        
        last_proj = *overlap_axis;
        
        // look in a direction normal to that face
        // and find a new point on the minkowski difference.
        Vec<T,N> minkowski_pt = shape_a.convexSupport( best_face->n) - 
                                shape_b.convexSupport(-best_face->n);
        
        // finish if we add a duplicate point (possible for polytopes).
        // such a point would produce a degenerate face, i.e. we aren't going
        // to find any more new faces. this is a linear algorithm-- but I am 
        // not sure whether a hashset construction/test is actually worth it.
        for (auto v : verts) if (v == minkowski_pt) goto FINISH;
        
        verts.push_back(minkowski_pt);
    
            // debug: emit newest hull pt
            EMIT(minkowski_pt);
        
        // now delete all the faces that fall "behind" the new point.
        for (auto f = faces.begin(); f != faces.end();) {
            // if this face is behind the new point...
            if (f->n.dot(minkowski_pt - verts[f->v[0]]) > 0) {
                // add all the edges of this face to the boundary edge list.
                // (all the non-internal edges are hole boundary, which must each
                //  become part of a new face).
                for (index_t i = 0; i < N; i++) {
                    // find a vertex to exclude from the face, making an edge.
                    detail::Edge<T,N> e;
                    for (index_t j = 0; j < N - 1; j++) {
                        e.v[j] = f->v[(i + j + 1) % N];
                    }
                    
                    // add the edge.
                    // if there exists an edge in the list with opposite winding, they "annihilate"
                    // but edges are only ever shared by exactly two faces, so we don't need to check the parity!
                    detail::Edge<T,N> e_tmp;
                    bool annihilated = false;
                    for (auto e_other = edges.begin(); e_other != edges.end(); ++e_other) {
                        e_tmp = *e_other;
                        if (std::abs(detail::simplex_comparison(e.v, e_tmp.v, N-1)) == 1) {
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
        
        for (auto e : edges) {
            EMIT(e, true);
        }
        
        // patch up the hole by constructing new faces,
        // connecting each hole-adjacent edge to the new vertex.
        index_t new_vert = verts.size() - 1;
        for (auto e = edges.begin(); e != edges.end(); ++e) {
            detail::Face<T,N> f;
            std::copy(e->v, e->v + N - 1, f.v);
            f.v[N-1] = new_vert;
            
            // compute face normal
            for (index_t i = 0; i < N - 1; i++) {
                splex_vtx_buf[i] = verts[f.v[i]] - minkowski_pt;
            }
            f.n = orthogonal(splex_vtx_buf).unit();
            
            // correct the face winding, if necessary.
            if (true or N > 3) {
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
    
    FINISH:
    
    if (emit_debug) std::cout << "=\n";
    
    return true;
}


} // end namespace geom

#endif	/* INTERSECT_H */

