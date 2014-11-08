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
        detail::StaticQueue<SS, detail::gjk_queue_size_fac(N) + 1> q;

        Vec<T,N> A = simplex->pts[simplex->n - 1];
        SS boundary;

        q.create_back()->createFrom(*simplex);

        // process the "frontier"
        while (q.size > 0) {
            SS &s = q.peek_front(); // note that `s` includes `A` implicitly.
            bool parent_is_boundary = true;
            for (index_t i = 0; i < s.n; i++) {
                SS &child = q.create_back()->createFrom(s, i, A);
                Vec<T,N> normal = child.getNormal();

                if (A.dot(normal) < 0) { 
                    // the origin is "outside" this sub-simplex, and therefore
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

} // namespace detail
 
// todo: template GJK entirely over the input shape type.
// todo: what if I want to do stuff like dilate or erode the minkowski difference?
//       what if I want to do stuff to only one of the shapes?
// todo: create intersect() function suite which defaults to GJK, but uses
//       cleverer things if they exist.
//       i.e. Intersector<A,B> uses GJK, but specializations may be written
//       for particular combinations of A and B.


// Use the Gilbert-Johnson-Keerthi algorithm to test for overlap between
// two convex shapes, and return the axis of overlap.
template <typename T, index_t N>
bool gjk_intersect(const Convex<T,N> &shape_a, 
                   const Convex<T,N> &shape_b,
                   Vec<T,N> *overlap_axis=NULL) {
    // choose an arbitrary initial direction.
    Vec<T,N> initial;
    if (overlap_axis and *overlap_axis != Vec<T,N>::zeros) 
        initial = *overlap_axis;
    else 
        initial[0] = 1;
    
    Vec<T,N> a = shape_a.convexSupport(initial) - 
                 shape_b.convexSupport(-initial);
    Vec<T,N> d = -a;
    detail::Simplex<T,N> s;
    s.insert(a);
    
    // supported by empirical measurement, plus a safe margin:
    const index_t cutoff = 10 * (1 << (N - 2)); 
    index_t i = 0;
    
    while (true) {
        a = shape_a.convexSupport( d) - 
            shape_b.convexSupport(-d);
        if (a.dot(d) < 0) return false;
        s.insert(a);
        
        if (detail::gjk_simplex_nearest_origin(&s, &d)) {
            if (overlap_axis) *overlap_axis = d;
            return true;
        }
        if (++i > cutoff) { return true; }
    }
}


template <typename T, index_t N>
inline bool gjk_intersect(const Vec<T,N> *pts_a, index_t n_a,
                          const Vec<T,N> *pts_b, index_t n_b,
                          Vec<T,N> *overlap_axis=NULL) {
    detail::UnstructuredPointcloud<T,N> a(pts_a, n_a);
    detail::UnstructuredPointcloud<T,N> b(pts_b, n_b);
    return gjk_intersect(a, b, overlap_axis);
}


} // end namespace geom

#endif	/* INTERSECT_H */

