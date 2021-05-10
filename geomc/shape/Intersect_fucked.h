/* 
 * File:   Intersect.h
 * Author: tbabb
 *
 * Created on October 14, 2014, 11:30 PM
 */

#ifndef INTERSECT_H
#define INTERSECT_H

#include <climits>
#include <geomc/linalg/Vec.h>
#include <geomc/linalg/Orthogonal.h>
#include <geomc/shape/Bounded.h>

#include <geomc/shape/Simplex.h>

#include <vector>
 
// xxx: this needs to be refactored to use mainstream Simplex<T,N>.
//      that implementation uses the stack instead of that god-awful 
//      StaticQueue garbage.
 
// todo: future: if we can raycast against the minkowski volume for two
//       concave shapes, then we can intersection test:
//       simply raycast away from the origin in both directions,
//       and count the parity of the intersections.
 
// xxx debug
#define EMIT(...) if(emit_debug) emit(__VA_ARGS__)

// todo: can we get faster/more stable/simpler results by solving for barycentric coords?
// todo: disjoint_separating_axis()
//       x occasional infinite loop in polytope explosion
//         > caused by final point being a duplicate, resulting in degenerate simplex
//       - ~25% of disjoint boxes get the wrong separation axis
//       - repro cases for both by running `distcheck` in `noodle`. see seeds below.
//       * must implement visualization protocol in disjoint_...()
//         - pay special scrutiny to backfacing case.
// todo: clean up all the debug crap
 
 
// alert!
// see: http://stackoverflow.com/questions/31738959/finding-the-point-on-a-simplex-closest-to-the-origin-using-gjks-distance-subalg
// there is an algorithm to find the closest point on a simplex in O(n) time.
// that will give some coordinate params for each of the bases, and then we can remove the bases that are zero.
// this is potentially much faster than the ad-hoc crap you are doing.


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
    
    
    // Store a sub-simplex, keeping track of its null basis and spanning basis.
    // The final vertex ("A") is implictly translated to the origin to simplify computation.
    template <typename T, index_t N>
    struct SubSimplex {
        
        SubSimplex():n(0) {}
        
        index_t pts[N+1]; // indecies of the verts into the source simplex
                          // in principle, this is not needed, as we could just 
                          // translate by `A`, but for floating point raisins
                          // we keep track of this.
                          // the final index is a point on the simplex.
                          // in intersection code, this is `A`.
        
        // basis for spanned space and its orthogonal complement. laid out like:
        //              parent         simplex
        //  normal    null basis    spanning basis
        //    [1]    [ N - n - 1 ]      [ n ]
        // All spanning bases point from `A` toward the corresponding vertex.
        Vec<T,N> basis[N];
        
        // number of stored points (i.e., excludes "A"). equal to the dimension spanned by `this`.
        index_t n;
        
        inline       Vec<T,N>* nullspace()               { return basis; }
        inline const Vec<T,N>* nullspace()         const { return basis; }
        inline       Vec<T,N>* spanning()                { return basis + N - n; }
        inline const Vec<T,N>* spanning()          const { return basis + N - n; }
        inline       Vec<T,N>& spanning(index_t i)       { return basis[N - n + i]; }
        inline const Vec<T,N>& spanning(index_t i) const { return basis[N - n + i]; }
        
        
        // create a new, complete simplex from `s`. The final vertex is of `s` is "A".
        SubSimplex& createFrom(const Simplex<T,N>& s) {
            n = s.n - 1;
            for (index_t i = 0; i < s.n; i++) pts[i] = i;
            
            // compute spanning and null spaces
            if (n > 0) {
                Vec<T,N> A = s[n];
                Vec<T,N> *spanning_basis = spanning();
                for (index_t i = 0; i < n; i++) {
                    spanning_basis[i] = s[i] - A;
                }
                // (short circuits if n == N)
                geom::nullspace(spanning_basis, n, basis);
            }
            return *this;
        }
        
        
        // create a new sub-simplex of lower dimension by excluding the vertex with index `excluded`.
        SubSimplex& createFrom(const SubSimplex& s, index_t excluded) {
            n = s.n - 1;
            
            // copy the source's null basis, leaving a hole for the normal to come
            std::copy(s.nullspace(), s.nullspace() + N - s.n, nullspace() + 1);
            
            // copy the source's spanning basis, skipping the excluded vtx.
            Vec<T,N>* spanning_basis = spanning();
            if (excluded < s.n) {
                index_t j = 0;
                for (index_t i = 0; i < s.n; i++) {
                    if (i == excluded) continue;
                    pts[j]            = s.pts[i];
                    spanning_basis[j] = s.spanning(i);
                    ++j;
                }
                pts[n] = s.pts[s.n]; // a pt on the simplex; use A.
            } else {
                // refactor the spanning basis to exclude "A"
                // (only used in contact-finding; ordinary EPA can, by construction,
                // skip checking sub-simplexes not involving "A").
                Vec<T,N> B = s.spanning(s.n - 1);
                for (index_t i = 0; i < n; i++) {
                    spanning_basis[i] = s.spanning(i) - B;
                    pts[i] = s.pts[i];
                }
                pts[n] = s.pts[n]; // pts[n] is a point on the simplex, used later.
                                   // we are exluding s.pts[s.n] (A), so use B instead.
            }
            // `basis` now looks like: [1 blank] [parent null space] [spanning basis]
            
            // compute the normal
            if (n > 0) {
                // the normal lies orthogonal to the null space of the parent
                // and orthogonal to the space of the child. the `n-1` tail of `basis` 
                // now contains both.
                Vec<T,N> normal = orthogonal(basis + 1);
                // flip the normal if it points "inside"
                Vec<T,N> v = (excluded == s.n) ? s.spanning(0) : s.spanning(excluded);
                if (normal.dot(v) > 0) normal = -normal;
                // extend the null basis.
                basis[0] = normal;
            } else {
                // we are a point; parent simplex is a line; there is only one normal
                basis[0] = -s.spanning(0);
            }
            return *this;
        }
        
        
        inline Vec<T,N> getNormal() const {
            return basis[0];
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
            if (size > N) throw "BAD";// xxx debug
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


// =============== begin debug {

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
void emit(const Vec<T,N>& v, bool highlight=false) {
    if (highlight) std::cout << "h ";
    std::cout << "v ";
    for (index_t i = 0; i < N; i++) std::cout << v[i] << " ";
    std::cout << "\n";
}

template <typename T, index_t N>
void emit(const detail::Simplex<T,N>& s, bool full) {
    if (full) {
        for (int i = 0; i < s.n; i++) {
            emit(s.pts[i]);
        }
    }
    detail::Edge<T,N> e;
    detail::Face<T,N> f;
    switch(s.n) {
        case 2: 
            // emit a line segment
            e.v[0] = 0;
            e.v[1] = 1;
            emit(e); 
            break;
        case 3:
            // emit a single face
            for (int i = 0; i < 3; i++) f.v[i] = i;
            emit(f);
            break;
        case 4:
            // emit 4 faces; excluding k each time
            for (int k = 0; k < 4; k++) {
                int ii = 0;
                for (int i = 0; i < 4; i++) {
                    if (i == k) continue;
                    f.v[ii++] = i;
                }
                emit(f);
            }
            break;
    }
}

// } =============== end debug
    
    
    constexpr index_t gjk_queue_size_fac(index_t n) {
        return (n < 2) ? 1 : n * gjk_queue_size_fac(n-1);
    }
    
    
    template <typename T, index_t N>
    struct SubSimplex_ {
        const detail::Simplex<T,N>& sup;
        index_t pts[N];
        index_t n;
        
        SubSimplex_(const detail::Simplex<T,N>& s, index_t exclude) : sup(s), n(s.n - 1) {
            index_t ct = 0;
            for (index_t i = 0; ct < n; ++i) {
                if (i == exclude) continue;
                pts[ct++] = i;
            }
        }
        
        SubSimplex_(const detail::SubSimplex_<T,N>& s, index_t exclude) : sup(s.sup), n(s.n - 1) {
            index_t ct = 0;
            for (index_t i = 0; ct < n; ++i) {
                if (i == exclude) continue;
                pts[ct++] = s.pts[i];
            }
        }
        
        inline const Vec<T,N>& operator[](index_t i) const {
            return sup[pts[i]];
        }
    };
    
    // never thought I'd have use for template templates.
    // this allows dj() to accept either a Simplex or a SubSimplex_. 
    // (construction of a sub-simplex from either is syntactically identical).
    // we can't just make it a single, named template type, because then we can't extract T and N.
    template <template <typename, index_t> class SimplexType, typename T, index_t N>
    T dj(const SimplexType<T,N>& s, index_t j) {
        if (s.n <= 1) return (T)1;
        detail::SubSimplex_<T,N> sj(s, j);
        
        T sum = (T)0;
        Vec<T,N> v = sj[0] - s[j];
        for (index_t i = 0; i < sj.n; ++i) {
            sum += dj(sj, i) * v.dot(sj[i]);
        }
        return sum;
    }
    
    
    // problems: 
    // - will often fail when a sub-simplex is selected (d will be wrong)
    // - rarely will fail after a repeat; usually picking an edge when a face should be picked.
    template <typename T, index_t N>
    Vec<T,N> gjk_johnson_simplex_nearest_origin(detail::Simplex<T,N>* s, bool* degenerate) {
        *degenerate = false;
        T sum = (T)0;
        detail::Simplex<T,N> out;
        Vec<T,N> p;
        for (index_t i = 0; i < s->n; ++i) {
            T d_j = dj(*s, i);
            if (d_j > 0) out[out.n++] = (*s)[i];
            p   += (*s)[i] * d_j;
            sum += d_j;
        }
        if (sum == 0) *degenerate = true;
        *s = out;
        return p / sum;
    }
    
    
    template <typename T, index_t N>
    bool gjk_direction_to_origin(detail::SubSimplex<T,N> &boundary, const Vec<T,N> &A, Vec<T,N> *d) {
        // find the direction to the origin. 
        // assumes a simplex which has already passed containment tests.
        // project the direction to the origin onto the null space of our simplex.
        // (the final `else` would be sufficient by itself, but in some cases
        // we can be clever and save some redundant calculations).
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
            // check for degenerate case: origin lies on face
            // (we do not check this in the other cases because `d` will represent a vector
            // from the origin to the simplex, and we check this against zero at the end).
            if (-A.projectOn(*d).mag2() == 0) return true;
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
                Vec<T,N> *proj_basis = boundary.nullspace();

                if ( null_dim > boundary.n ) {
                    proj_basis = boundary.spanning();
                    proj_n = boundary.n;
                    swap_proj = true;
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
                Vec<T,N> vs[3] = {(T)0, -boundary.spanning(0), -A};
                vs[0] = orthogonal(vs+1);
                *d = orthogonal(vs);
            }
        }
        // if the distance from the simplex to the origin is zero,
        // then the simplex "contains" the origin (i.e. the origin lies on it).
        // we must terminate now, because further simplexes will be degenerate.
        return d->mag2() == 0;
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
     *
     * edit: are we simply guaranteed to terminate before we check more than N edges?
     */
    template <typename T, index_t N>
    bool gjk_simplex_nearest_origin(detail::Simplex<T,N> *simplex, Vec<T,N> *d, bool fullcheck=false) {
        typedef detail::SubSimplex<T,N> SS;
        // queue for sub-simplexes to check.
        // (sub-simplexes store their spanning and null bases).
        // xxx: This is questionable. replace with deque / abstract to cache object.
        //      ...or can we factor to recursion, and use the stack?
        detail::StaticQueue<SS, 2 * detail::gjk_queue_size_fac(N) + 1> q; // xxx doubled in size :S
        
        SS boundary;
        
        // enqueue this simplex, and its null basis 
        q.create_back()->createFrom(*simplex);

        // process the "frontier"
        while (q.size > 0) {
            SS &parent = q.peek_front();
            // determine whether the origin projects onto this simplex, or instead to one of its edges.
            // edge simplexes must be recursively checked for projection containment and so are enqueued.
            bool parent_is_boundary = true;
            index_t ct = (fullcheck and parent.n > 0) ? parent.n + 1 : parent.n;
            for (index_t i = 0; i < ct; i++) {
                SS &child = q.create_back()->createFrom(parent, i);
                Vec<T,N> normal = child.getNormal();
                Vec<T,N> B = simplex->pts[child.pts[0]]; // a pt on the simplex
                
                if (B.dot(normal) < 0) { 
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
                boundary = parent;
                break;
            } else {
                // this simplex is checked; pop it off the stack.
                q.destroy_front();
            }
        }
        
        bool containment = gjk_direction_to_origin(boundary, simplex->pts[boundary.pts[boundary.n]], d);
        
        // put the new simplex back into the return variable.
        // (boundary's indecies are strictly increasing, so this is ok).
        for (index_t i = 0; i < boundary.n + 1; i++) {
            (*simplex)[i] = (*simplex)[boundary.pts[i]];
        }
        simplex->n = boundary.n + 1;
        
        return containment;
    }
    
    template <typename T, index_t N>
    struct UnstructuredPointcloud : virtual public Convex<T,N> {
        const Vec<T,N>* pts;
        index_t n;
        
        UnstructuredPointcloud(const Vec<T,N> *pts, index_t n):pts(pts),n(n) {}
        
        Vec<T,N> convex_support(Vec<T,N> d) const {
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
                       detail::Simplex<T,N> *s,
                       bool emit_debug=false) {
        // choose an arbitrary initial direction.
        Vec<T,N> initial;
        if (overlap_axis and *overlap_axis != Vec<T,N>::zeros) 
            initial = *overlap_axis;
        else 
            initial[0] = 1;
        
        Vec<T,N> a = shape_a.convex_support( initial) - 
                     shape_b.convex_support(-initial);
        Vec<T,N> d = -a;
        s->n = 0;
        s->insert(a);
        
        // supported by empirical measurement, plus a safe margin:
        const index_t cutoff = 10 * (1 << (N - 2)); 
        index_t i = 0;
        
        while (true) {
            a = shape_a.convex_support( d) - 
                shape_b.convex_support(-d);
            
            s->insert(a);
            
            if (a.dot(d) < 0) {
                if (overlap_axis) *overlap_axis = -a.projectOn(d);
                return false;
            }
            
            if (detail::gjk_simplex_nearest_origin(s, &d)) {
                // the "nearest" simplex containing a projection of the origin is a complete volume, 
                // or the origin lies on its projection
                if (overlap_axis) *overlap_axis = d;
                return true;
            }
            
            if (++i > cutoff) { return true; }
        }
    }
    
    
    template <typename T, index_t N>
    bool gjk_intersect_johnson(const Convex<T,N> &shape_a, 
                       const Convex<T,N> &shape_b,
                       Vec<T,N> *overlap_axis,
                       detail::Simplex<T,N> *s,
                       bool emit_debug=false) {
        // choose an arbitrary initial direction.
        Vec<T,N> initial;
        if (overlap_axis and *overlap_axis != Vec<T,N>::zeros) 
            initial = *overlap_axis;
        else 
            initial[0] = 1;
        
        Vec<T,N> a = shape_a.convex_support( initial) - 
                     shape_b.convex_support(-initial);
        Vec<T,N> d = -a;
        s->n = 0;
        s->insert(a);
        
        // supported by empirical measurement, plus a safe margin:
        const index_t cutoff = 10 * (1 << (N - 2)); 
        index_t i = 0;
        
        while (true) {
            a = shape_a.convex_support( d) - 
                shape_b.convex_support(-d);
            
            s->insert(a);
            
            // xxx debug {
            if (emit_debug) {
                std::cout << s->n << " ";
                for (index_t k = 0; k < s->n; ++k) {
                    std::cout << (*s)[k];
                }
                std::cout << a << "\n";
            }
            // }
            
            if (a.dot(d) < 0) {
                if (overlap_axis) *overlap_axis = -a.projectOn(d);
                return false;
            }
            
            index_t old_n = s->n;
            bool degenerate = false;
            Vec<T,N> new_d = -detail::gjk_johnson_simplex_nearest_origin(s, &degenerate);
            if (s->n < old_n) {
                // xxx hack
                // for some reason things fail when simplex is a subset.
                new_d = -detail::gjk_johnson_simplex_nearest_origin(s, &degenerate);
            }
            
            // xxx debug {
            if (emit_debug) {
                std::cout << s->n << " ";
                for (index_t k = 0; k < s->n; ++k) {
                    std::cout << (*s)[k];
                }
                std::cout << -new_d << "\n";
            }
            // }
            
            if (degenerate) {
                // we failed to "escape" the previous face
                // there cannot be a point closer to the origin; there is no intersection.
                if (overlap_axis) *overlap_axis = d;
                return false;
            }
            
            if (s->n > N or new_d == Vec<T,N>::zeros) {
            // if (detail::gjk_simplex_nearest_origin(s, &d)) {
                // xxx handle export of `d` in case of d = 0
                // the "nearest" simplex containing a projection of the origin is a complete volume, 
                // or the origin lies on its projection
                if (overlap_axis) *overlap_axis = d;
                return true;
            }
            
            d = new_d;
            if (++i > cutoff) { return true; }
        }
    }
    
    
    template <typename T, index_t N>
    void explode_simplex(const Convex<T,N>& shape_a,
                         const Convex<T,N>& shape_b, 
                         detail::Simplex<T,N>* splex) {
       while (splex->n < N + 1) {
            // search the initial nullspace axes first, and only bother to
            // re-construct the sub-simplex and its nullspace if those searches failed.
            detail::SubSimplex<T,N> s_initial;
            s_initial.createFrom(*splex);
            index_t n_null = N - s_initial.n;
            // search in both directions
            for (index_t i = 0; i < n_null * 2 and splex->n < N + 1; i++) {
                index_t ii = i % n_null;
                T      sgn = (i / n_null) ? ((T)-1) : ((T)1);
                Vec<T,N> n = sgn * s_initial.nullspace()[ii];
                Vec<T,N> a = shape_a.convex_support( n) - 
                             shape_b.convex_support(-n);
                // do not insert a duplicate vtx.
                bool dupe = false;
                for (index_t j = 0; j < splex->n; j++) {
                    if ((*splex)[j] == a) {  // todo: O(n) search ok?
                        dupe = true;
                        break;
                    }
                }
                if (not dupe) {
                    splex->insert(a);
                }
            }
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
                   Vec<T,N> *overlap_axis=NULL,
                   bool emit_debug=false) {  // xxx debug
    detail::Simplex<T,N> s;
    return detail::gjk_intersect(shape_a, shape_b, overlap_axis, &s);
}


template <typename T, index_t N>
inline bool gjk_intersect(const Vec<T,N> *pts_a, index_t n_a,
                          const Vec<T,N> *pts_b, index_t n_b,
                          Vec<T,N> *overlap_axis=NULL,
                          bool emit_debug=false) { // xxx debug
    detail::UnstructuredPointcloud<T,N> a(pts_a, n_a);
    detail::UnstructuredPointcloud<T,N> b(pts_b, n_b);
    return gjk_intersect(a, b, overlap_axis, emit_debug);
}


// todo: write an optimized N=2 implementation and then verify it against the general one
//       - edge comparison is a hard equality check
//       - there can only ever be exactly two unpatched hole edges.
// todo: optimize face insertion
// todo: abstract into class to allow buffer amortization
// todo: test at what point hashtable construction is faster for duplicate pt tests.
// todo: keep the faces in sorted order by distance. insertion is O(lg(n)) and
//       searching is O(1). std::list means no shuffling.
// todo: abstract penetration solving into a class, with cache workspace
//       holding the std::lists, etc.
// todo: clean up debug mumbo jumbo

// observation: The number of faces is strictly increasing.
//              you don't need to delete a face [O(n)]; they can be clobbered.
//              use a ~~skip list~~ ~~std::list~~ sorted multimap.

// xxx: inf loop case is `distcheck 18292`

// xxx: disjoint shapes do not find the correct separation axis ~25% of the time.
//    - the problem is:
//      - instead of adding A to the nearest simplex (which may be a pt, e.g.)
//        you make a new volume simplex using A.
//        - assuming this is sufficient, rather than maintaining a full boundary manifold
//      x simplex_nearest_origin() picks the *first* simplex which
//        completely "contains" the origin. for intersection code, this is sufficient,
//        because there can only be one new fully-containing simplex (if they are checked
//        in hierarchical order). when we begin from a full tetrahedron this is not so
//        because the extrusions of three faces can make a pyramid on the "backside"
//        of the simplex. we need a full distance-based check, and formulating this to be
//        O(1) is a tricky endeavor.
//        > no, this isn't it.
//
//    - (seed) `distcheck 18299` repros this.
//    x apparently backface checking was just never turned on.
//      > this is fixed
//    x during simplex_nearest_origin we eventually construct a simplex of size zero, then
//      try to make a sub-simplex of lower dimension from it.
//      > this is fixed
//
//    - d is zero in some failure cases.
//    - the problem manifests as gjk_direction_to_origin erroneously reporting full volume containment of the origin.
//
// 
//    old problems, from before we enabled backface checking:
//    - there may be an issue in that the notion of "backfacing" is destroyed when "A" is discarded on recursion?
//    - consider the sign of the normal <<< 
//    - consider projecting origin to simplex null basis:
//      - concerned about case where said projection is zero, e.g.:
//        an edge points directly at the origin, as when the minkowski volume is a sphere. 

// xxx: profile test of average separation error is wrong-o.
//      convex_support(d) does not give you the contact point.
//      you need the individual support points from the two shapes.

// xxx: new inf loop in disjoint test.

// given a minkowski volume which does *not* enclose the origin, find the closest point
// on the volume to the origin.
template <typename T, index_t N>
void disjoint_separation_axis(const Convex<T,N>& shape_a,
                              const Convex<T,N>& shape_b, 
                              Vec<T,N> *overlap_axis,
                              detail::Simplex<T,N>* splex,
                              double fractional_tolerance = 0.001,
                              index_t iteration_limit = -1,
                              bool emit_debug = false) {
    
    // final pt may have been a duplicate; if so, remove it.
    Vec<T,N> A = splex->pts[splex->n - 1];
    for (index_t i = 0; i < splex->n - 1; i++) {
        if (splex->pts[i] == A) {
            splex->remove(splex->n - 1);
            break;
        }
    }
    
    // expand the simplex to a volume
    detail::explode_simplex(shape_a, shape_b, splex);
    
    // iteratively search outward from the closest (sub-)simplex.
    Vec<T,N> last_proj;
    Vec<T,N> new_proj;
    bool looped = false;
    for (index_t j = 0; j != iteration_limit; j++) {
        Vec<T,N> d;
        
        // update `splex` and `d`
        gjk_simplex_nearest_origin(splex, &d, true);
        
        // search in the direction of the "best" normal
        Vec<T,N> a = shape_a.convex_support( d) -
                     shape_b.convex_support(-d);
        
        // project the origin onto the simplex, down along its normal
        new_proj = ((*splex)[0]).projectOn(d);
        
        // if the search point already exists in the simplex, halt
        for (index_t i = 0; i < splex->n; i++) if ((*splex)[i] == a) {
            goto DONE_DJSA;
        }
    
        // if our estimated closest point is not sufficiently different from
        // our last estimate (or if the origin is on the current face),
        // decide that we've converged and quit.
        if ((looped and fractional_tolerance > 0 and 
                (last_proj - new_proj).mag() / new_proj.mag() < fractional_tolerance) 
                // origin is exactly on the face:
                or new_proj == Vec<T,N>::zeros) {
            break;
        }
        
        splex->insert(a);
        last_proj = new_proj;
        looped = true;
    }
    DONE_DJSA: *overlap_axis = new_proj;
}


/**
 * @ingroup shape
 * Use the Expanding Polytope Algorithm to find a minimum translation vector that would bring shape B into contact with A.
 * If A and B are not interpenetrating the separation axis is undefined. 
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
        // xxx: this shit below ain't work
        disjoint_separation_axis(shape_a, shape_b, 
                                 overlap_axis,
                                 &splex,
                                 fractional_tolerance,
                                 iteration_limit);
        return false;
    }
    
    if (splex.n != N + 1) {
        // rare degenerate case: origin lies on a sub-simplex
        // simplex must be "exploded" to a volume.
        detail::explode_simplex(shape_a, shape_b, &splex);
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
        
        // find the face closest to the origin
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
        Vec<T,N> minkowski_pt = shape_a.convex_support( best_face->n) - 
                                shape_b.convex_support(-best_face->n);
        
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
        
        // debug
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
    
    FINISH:
    
    if (emit_debug) std::cout << "=\n";
    
    return true;
}


} // end namespace geom

#endif  /* INTERSECT_H */

