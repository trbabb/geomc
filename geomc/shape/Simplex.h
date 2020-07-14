/* 
 * File:   Frustum.h
 * Author: tbabb
 *
 * Created on November 8, 2014, 11:09 PM
 */

#ifndef SIMPLEX_H
#define	SIMPLEX_H

#include <geomc/shape/Bounded.h>
#include <geomc/linalg/Vec.h>
#include <geomc/linalg/Orthogonal.h>

// todo: orthogonal projection to a subspace is a very general linalg operation,
//       and projection_contains() makes use of it. The operation that does this is
//       the pseudoinverse, which can be computed efficiently using a QR decomposition.
//       Probably the best way to do QR is with Householder transformations; the algorithm
//       can be found complete here: http://people.inf.ethz.ch/gander/papers/qrneu.pdf
//       
//       the above should probably be exposed as orthogonal_project(basis[], n, P) in
//       Orthogonal.h, and invoked here by projection_contains(). Also consider adding
//       subspace_projection(bases[], n) -> Mtx, which can amortize the QR bidness.
//
//       right now we are running two LUP decompositions, and this is likely
//       inefficient/unstable.

// todo: handle degeneracy. GJK can make use of it.
// todo: performance check old GJK vs new
// todo: performance check projection_contains() vs. project() == this


namespace geom {

/**
 * @ingroup shape
 * @brief A simplex in N dimensions (e.g. a tetrahedron, triangle, line, point).
 *
 * The simplex may contain up to `N + 1` points, in which case it encloses a 
 * volume, if it is not degenerate. If it contains fewer points, then it
 * spans a subspace of the space in which it is embedded.
 */
template <typename T, index_t N>
class Simplex : public virtual Convex<T,N> {

public:

    /// Vertices of this simplex.
    Vec<T,N> pts[N + 1];
    /// Number of vertices in this simplex.
    index_t n;
    

    /// Construct an empty simplex, with no vertices.
    Simplex():n(0) {}

    /// Construct an `n`-cornered simplex with vertices at `verts`.
    Simplex(const Vec<T,N>* verts, index_t n):n(n) {
        std::copy(verts, verts + n, pts);
    }


    /// Get the `i`th vertex in this simplex.
    Vec<T,N>& operator[](index_t i) {
        return pts[i];
    }


    /// Get the `i`th vertex in this simplex.
    Vec<T,N> operator[](index_t i) const {
        return pts[i];
    }
    
    
    /// Simplex equality check.
    bool operator==(const Simplex<T,N>& other) const {
        // overridden because there may be garbage in the 
        // unused vertices which shouldn't count as a difference.
        if (other.n != n) return false;
        for (index_t i = 0; i < n; ++i) {
            if (other.pts[i] != pts[i]) return false;
        }
        return true;
    }
    
    
    /// Simplex inequality check.
    inline bool operator!=(const Simplex<T,N>& other) const {
        return not ((*this) == other);
    }


    /**
     * @brief Return a copy of this simplex with an additional vertex at `p`.
     *
     * If this simplex already has `N + 1` points, a copy of this simplex is returned.
     */
    Simplex<T,N> operator|(const Vec<T,N>& p) const {
        Simplex<T,N> s = *this;
        if (n < N + 1) {
            s.n += 1;
            s.pts[n] = p;
        }
        return s;
    }


    /**
     * @brief Extend this simplex to include `p` by adding `p` as a vertex.
     *
     * If this Simplex already has `N+1` vertices, then this operator has no effect.
     */
    Simplex<T,N>& operator|=(const Vec<T,N>& p) {
        if (n < N + 1) {
            pts[n] = p;
            n += 1;
        }
        return *this;
    }
    
    /**
     * @brief Apply a transformation to the points of this Simplex.
     */
    Simplex<T,N>& operator*=(const AffineTransform<T,N>& xf) {
        for (index_t i = 0; i < n; ++i) {
            pts[i] = xf * pts[i];
        }
        return *this;
    }
    
    /**
     * @brief Apply an inverse transformation to the points of this Simplex.
     */
    Simplex<T,N>& operator/=(const AffineTransform<T,N>& xf) {
        for (index_t i = 0; i < n; ++i) {
            pts[i] = pts[i] / xf;
        }
        return *this;
    }
    
    /**
     * @brief Project the point to the simplex along its orthogonal basis
     * (if there is one) and test for containment.
     * 
     * @param p A point.
     * @param params Optional return variable; return the `n - 1` surface parameters
     * of this simplex representing the projection of `p` onto the space
     * spanned by this simplex.
     * @return `true` if `p` is on or inside this simplex after projection; 
     * `false` otherwise.
     */
    bool projection_contains(const Vec<T,N>& p, Vec<T,N>* params=nullptr) const {
        switch (n) {
            case 0:
                return false;
            
            case 1: {
                return true;
            } break;
            
            case 2: {
                Vec<T,N> v = pts[1] - pts[0];
                T d = v.dot(p - pts[0]) / v.mag();
                if (params) (*params)[0] = d;
                return d >= 0 and d <= 1;
            } break;
            
            default: {
                // simplex spans a subspace. orthogonally project p along its nullspace.
                // if the simplex spans the domain, then nullspace() early exits and
                // everything else works as normal.
                
                // nullspace comes first (cause we won't solve for its coords),
                // followed by the simplex basis.
                const index_t n_bases = n - 1;
                const index_t n_null  = N - n + 1;
                Vec<T,N> null_bases[N];
                Vec<T,N>* const bases = null_bases + n_null; 
                
                // compute simplex bases
                for (index_t i = 1; i < n; ++i) {
                    bases[i - 1] = pts[i] - pts[0];
                }
                
                // find the simplex's nullspace
                nullspace(bases, n_bases, null_bases);
                
                // write P as a linear combination of the simplex basis and 
                // its nullspace bases. (but don't bother actually solving for 
                // the nullpsace coords; we don't use them).
                Vec<T,N> x;
                if (not linear_solve(bases, &x, p - pts[0], n_null)) return false;
                
                // check that the coordinates of `p` are inside the simplex
                T sum = 0;
                for (index_t i = n_null; i < n_bases; ++i) {
                    if (x[i] < 0 or sum > 1) return false;
                    sum += x[i];
                }
                
                // provide surface parameters to caller, e.g. so they can perform
                // the projection if they want to.
                if (params) {
                    for (index_t i = 0; i < n_bases; ++i) {
                        (*params)[i] = x[i + n_null];
                    }
                }
                
                return (sum <= 1);
            }
        }

        // statically unreachable, but some less-clever compilers may complain without:
        return false;
    }
    
    /**
     * @brief Simplex-point intersection test.
     * 
     * If this simplex has fewer than `N + 1` vertices, this simplex
     * is not a volume, and this function returns `false`. 
     * 
     * @param p A point.
     * @return `true` if `p` is on or inside this simplex; `false` otherwise.
     */
    bool contains(const Vec<T,N>& p) const {
        if (n < N + 1) return false;
        
        // solve for barycentric coordinates x:
        // Σ (1, v[i]) * x[i] = (1, p)
        
        const index_t K = N + 1;
        T m[K * K];
        T x[K];
        T b[K];
        for (index_t i = 0; i < K; ++i) {
            const T* v = pts[i].begin();
            m[i * K] = 1;
            std::copy(v, v + N, m + i * K + 1);
        }
        b[0] = 1;
        std::copy(p.begin(), p.end(), b + 1);
        
        linear_solve<T,false>(m, K, x, b, 1);
        // ^^ skip solving for the sum of the weights; we know it's 1.
        
        T sum = 0;
        for (index_t i = 1; i < K; ++i) {
            sum += x[i];
            if (x[i] < 0 or sum > 1) return false;
        }
        return true;
    }


    /**
     * @brief Create a new sub-simplex by excluding the `i`th vertex in this simplex.
     */
    Simplex<T,N> exclude(index_t i) const {
        Simplex<T,N> s;
        index_t c = 0;
        for (index_t j = 0; j < n; ++j) {
            if (i != j) s.pts[c++] = pts[j];
        }
        s.n = n - 1;
        return s;
    }

    
    /**
     * @brief Project `p` to the nearest point on the simplex.
     *
     * The point is unchanged if `p` is already inside the simplex.
     * 
     * @param p The point to project onto the surface of this simplex.
     * @param onto Optional output parameter to receive the sub-simplex
     * onto which `p` was projected. It is permissible for `onto` to 
     * alias `this`.
     */
    Vec<T,N> project(const Vec<T,N>& p, Simplex<T,N>* onto=nullptr) const {
        Vec<T,N> buffer[N];
        Simplex<T,N> s;
        
        // why the recursion? why not project as in projection_contains(),
        // obtaining surface parameters, and then clip them?
        // A: Because in general, the subspace of the simplex will not
        //    be orthogonal, so clipping the coordinates to the unit triangle
        //    will project the point non-orthogonally along the bases of 
        //    that subspace. Instead we have to first find the sub-simplex 
        //    which is closest, and then project orthogonally down its nullspace.
        
        // initialize the nullspace, if necessary
        if (n - 1 < N) {
            Vec<T,N>* const basis      = buffer + N - n + 1;
            Vec<T,N>* const null_basis = buffer;
            for (index_t i = 1; i < n; ++i) {
                basis[i - 1] = pts[i] - pts[0];
            }
            nullspace(basis, n - 1, null_basis);
        }
        
        // this puts the simplex which faces `p` into `s`, and
        // leaves `buffer` containing its nullspace and spanning space.
        this->_find_nearest(p, buffer, false, Vec<T,N>::zeros, &s);
        
        // choose the smaller basis to project on
        index_t   n_bases = s.n - 1;
        index_t   n_null  = N - n_bases;
        index_t   proj_dims;
        Vec<T,N>* proj_basis;
        Vec<T,N>  out;
        bool proj_to_null;
        if (n_bases < n_null) {
            proj_dims    = n_bases;
            proj_basis   = buffer + n_null;
            proj_to_null = true;
        } else {
            proj_dims    = n_null;
            proj_basis   = buffer;
            proj_to_null = false;
        }
        
        // project to the relevant subspace
        if (proj_dims > 0) {
            if (N > 3 and (not proj_to_null or N - n > 0) and proj_dims > 1) {
                // if we're projecting to the nullspace, it's already orthogonal,
                // unless the nullspace started with > 1 dimension (in which case 
                // it was made by `nullspace` above; orthogonalize it now). 
                // this will never happen in 3D because there is always a projection
                // space of dim <= 1 to choose.
                orthogonalize(proj_basis, proj_dims);
            }
            Vec<T,N> p_in_basis = p - s.pts[0];
            for (index_t i = 0; i < proj_dims; ++i) {
                out += p_in_basis.projectOn(proj_basis[i]);
            }
        }
        
        // offset the subspace projection so it lies on the simplex surface
        if (proj_to_null) {
            // `out` is the direction away from the simplex toward `p`.
            // travel backwards in that direction from `p` to get to the simplex.
            out = p - out;
        } else {
            // `out` is the direction from s's corner to the projected point.
            // add s's corner to get the absolute point.
            out += s.pts[0];
        }
        
        if (onto) *onto = s;
        
        return out;
    }
    
    
    Rect<T,N> bounds() const {
        Rect<T,N> bb;
        for (index_t i = 0; i < n; ++i) {
            bb |= pts[i];
        }
        return bb;
    }
    
    
    Vec<T,N> convex_support(Vec<T,N> d) const {
        T best    = pts[0].dot(d);
        index_t k = 0;
        for (index_t i = 1; i < n; ++i) {
            T a = pts[i].dot(d);
            if (a > best) {
                best = a;
                k = i;
            }
        }
        return pts[k];
    }

private:

    // `all_bases` shall be laid out like:
    //   [normal] [parent null basis] [spanning basis]
    // where `normal` remains to be computed.
    // in the base case (not is_sub), all nullspace axes are provided.
    bool _find_nearest(
            Vec<T,N> p, 
            Vec<T,N> all_bases[], 
            bool is_sub,
            Vec<T,N> missing_basis,
            Simplex<T,N>* onto) const {
        if (n == 1) {
            // projection onto a single point is trivial
            // note: we leave `all_bases` in a garbage state,
            // but we don't use it, since we don't need it!
            *onto = *this;
            return true;
        } else {
            // orthogonally project p onto this simplex along its nullspace.
            const index_t n_null       = N - n + 1;
            Vec<T,N>* const bases      = all_bases + n_null;
            Vec<T,N>* const null_bases = all_bases + 1;
            
            // compute simplex's basis
            for (index_t i = 1; i < n; ++i) {
                bases[i - 1] = pts[i] - pts[0];
            }
            
            // complete the simplex's nullspace
            Vec<T,N> normal;
            if (n_null != 0) {
                if (is_sub) {
                    // newest nullspace axis not provided; must be computed.
                    all_bases[0] = orthogonal(null_bases); // <- also contains `bases`! 
                
                    // `normal` is the vector pointing to `p` along the new axis
                    // (which will be orthogonal both to this simplex, and
                    // to the nullspace of the parent simplex).
                    normal = (p - pts[0]).projectOn(all_bases[0]);
                    
                    // might `p` be on the "inside" of the parent simplex?
                    if (normal.dot(missing_basis) > 0) {
                        // `p` falls "inward" of the boundary space spanned
                        // by this simplex. `p` therefore does not project to this
                        // simplex; instead it projects to some other simplex,
                        // possibly a parent of this, which *does* face (or envelop) `p`.
                        return false;
                    }
                }
                // since there is at least one nullspace axis,
                // move it aside for the sub-simplex searching code coming below,
                // which will prepend its own normal.
                bases[0] = all_bases[0];
            }
            
            // does `p` project to some sub-simplex of ours?
            // try each sub-simplex.
            for (index_t i = 0; i < n; ++i) {
                Simplex<T,N> sub = this->exclude(i);
                if (sub._find_nearest(p, all_bases, true, pts[i] - sub[0], onto)) {
                    // some sub-simplex of us both faces `p` and contains
                    // its projection. `all_bases` and `onto` now reflect that sub-simplex.
                    return true;
                }
            }
            
            // restore the state of `all_bases`
            // nobody messed with the normal because it was in their null space.
            // and the simplex which just filled in its basis was the one that 
            // excluded vertex `n`.
            all_bases[0] = bases[0];
            bases[0]     = pts[1] - pts[n];
            
            // no sub-simplex claimed `p` was "beyond" it and took the projection.
            // therefore, this is the simplex onto which `p` projects.
            *onto = *this;
            return true;
        }
    }


}; // class simplex



} // namespace geom

#endif	/* SIMPLEX_H */

