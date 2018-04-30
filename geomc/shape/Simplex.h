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
#include <geomc/linalg/LUDecomp.h>

// todo: Store A + basis vectors;
//       implement operator[] to do the subtraction.
//       translation would be cheap this way too (move only A!)

namespace geom {

/**
 * @ingroup shape
 * @brief A simplex in N dimensions (e.g. a tetrahedron, triangle, line, point).
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


        /**
         * @brief Return a copy of this simplex with an additional vertex at `p`.
         *
         * If this simplex already has `N + 1` points, a copy of this simplex is returned.
         */
        Simplex<T,N> operator+(const Vec<T,N>& p) const {
            Simplex<T,N> s = *this;
            if (n < N + 1) {
                s.n += 1;
                s[n] = p;
            }
            return s;
        }


        /**
         * @brief Add a new vertex to this simplex at `p`.
         */
        Simplex<T,N>& operator+(const Vec<T,N>& p) {
            if (n < N + 1) {
                pts[n] = p;
                n += 1;
            }
            return *this;
        }


        /**
         * @brief Project the point to the simplex along its orthoginal basis
         * (if there is one) and test for containment.
         * 
         * @param p A point.
         * @param params Optional return variable; return the `n - 1` surface parameters
         * of this simplex representing the projection of `p` onto the space
         * spanned by this simplex. (Variables in `params` with index of `n - 1`
         * or greater will contain nonsense values).
         * @return `true` if `p` is on or inside this simplex after projection; 
         * `false` otherwise.
         */
        bool projection_contains(const Vec<T,N>& p, Vec<T,N>* params=nullptr) const {
            switch (n) {
                case 0:
                    return false;
                
                case 1: {
                    return pts[0] == p;
                } break;
                
                case 2: {
                    Vec<T,N> v = pts[1] - pts[0];
                    T d = v.dot(p - pts[0]) / v.mag();
                    if (params) (*params)[0] = d;
                    return d > 0 and d < 1;
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
                    if (not linearSolve(bases, &x, p - pts[0], n_null)) return false;
                    
                    // provide surface parameters to caller, e.g. so they can perform
                    // the projection if they want to.
                    if (params) {
                        *params = x;
                    }
                    
                    // check that the coordinates of `p` are inside the simplex
                    T sum = 0;
                    for (index_t i = 0; i < n_bases; ++i) {
                        if (x[i] < 0) return false;
                        sum += x[i];
                    }
                    return (sum <= 1);
                }
            }

            // for pickier compilers:
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
            return (n - 1 == N) ? projection_contains(p) : false;
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
         * @param p The point to project onto the surface of this simplex.
         * @param onto Optional output parameter to receive the sub-simplex
         * onto which `p` was projected.
         *
         * The point is unchanged if `p` is already inside the simplex.
         */
        Vec<T,N> project(const Vec<T,N>& p, Simplex<T,N>* onto=nullptr) const {
            Vec<T,N> buffer[N];
            Vec<T,N> out;
            if (n < N + 1) {
                // we need to compute the nullspace.
                Vec<T,N>* const basis      = buffer + N - n + 1;
                Vec<T,N>* const null_basis = buffer;
                for (index_t i = 1; i < n; ++i) {
                    basis[i - 1] = pts[i] - pts[0];
                }
                nullspace(basis, n - 1, null_basis);
            }
            this->_project(p, buffer, false, out, &out, onto); // first `out` is dummy
            
            // todo: the buffer has your bases + nullspace in a good state.
            // consider adding option to solve for surface parameters.
            
            return out;
        }
        
        
        Vec<T,N> convexSupport(Vec<T,N> d) const {
            T best    = pts[0].dot(d);
            index_t k = 0;
            for (index_t i = 1; i < n; ++i) {
                T a = pts[i].dot(d);
                if (a > best) k = i;
            }
            return pts[k];
        }

    private:

        // `all_bases` shall be laid out like:
        //   [normal] [parent null basis] [spanning basis]
        // where `normal` remains to be computed.
        // in the base case (not is_sub), all nullspace axes are provided.
        bool _project(
                Vec<T,N> p, 
                Vec<T,N> all_bases[], 
                bool is_sub,
                Vec<T,N> missing_basis, 
                Vec<T,N>* proj, 
                Simplex<T,N>* onto) const {
            switch (n) {
                case 0:
                    // cannot project onto empty set
                    return false;
                
                case 1: {
                    // projection onto a single point is trivial
                    *proj = pts[0];
                    if (onto) *onto = *this;
                    return true;
                } break;
                
                default: {
                    // orthogonally project p onto this simplex along its nullspace.
                    const index_t n_bases = n - 1;
                    const index_t n_null  = N - n + 1;
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
                        
                            // the vector pointing to `p` along the new axis
                            // (which will be orthogonal both to this simplex, and
                            // to the nullspace of the parent simplex).
                            normal = (p - pts[0]).projectOn(all_bases[0]);
                            
                            // might `p` be on the "inside" of the parent simplex?
                            if (is_sub and normal.dot(missing_basis) > 0) {
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
                        if (sub._project(p, all_bases, true, pts[i] - sub[0], proj, onto)) {
                            // some sub-simplex of us both faces `p` and contains
                            // its projection. `proj` and `onto` now reflect that projection.
                            return true;
                        }
                    }
                    
                    // no sub-simplex claimed `p` was "beyond" it and took the projection.
                    // therefore, this is the simplex onto which `p` projects.
                    *proj = p - normal;
                    if (onto) *onto = *this;
                    return true;
                }
            }
        }


}; // class simplex
    
} // namespace geom

#endif	/* SIMPLEX_H */

