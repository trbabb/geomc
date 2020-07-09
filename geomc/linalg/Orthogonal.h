/* 
 * File:   Orthogonal.h
 * Author: tbabb
 *
 * Created on October 15, 2014, 12:23 AM
 */

#ifndef ORTHOGONAL_H
#define	ORTHOGONAL_H

#include <algorithm>
#include <geomc/linalg/Vec.h>
#include <geomc/linalg/LUDecomp.h>

/*
http://en.wikibooks.org/wiki/Linear_Algebra/Null_Spaces

get to row-echelon form, and solve:
 d e f g   x   0
 0 h i j * y = 0 
 0 0 k l   z   0
           w   1

       ↓   ↓
 d e f # g #   x   0
 0 h i # j # * y = 0
 0 0 k # l #   z   0
               -   0  ←┬─ ignore/clear extra cols/vars
               w   1   │
               -   0  ←┘

dx + ey + fz + g = 0
0x + hy + iz + j = 0
0x + 0y + kz + l = 0
0x + 0y + 0z + w = 1
--
          kz = -l
           z = -l / k
           ...
           etc.
 */

#define _MxElem(m,r,c) m[N*r + c]

namespace geom {
    
    // todo: orthogonal() does not respect winding.
    // the wedge product is maybe something you need.
    
    /**
     * @addtogroup linalg
     * @{
     */
    
    /**
     * Return a vector orthogonal to the given `N-1` vectors.
     * 
     * If any of the given basis vectors are linearly dependent,
     * the function returns the 0 vector.
     * 
     * @param v Array of `N-1` basis vectors.
     * @return A vector normal to all the members of `v`.
     */
    template <typename T, index_t N>
    Vec<T,N> orthogonal(const Vec<T,N> v[N-1]) {
        Vec<T,N> o;
        o[N-1] = 1;
        index_t P[N];
        T       m[N * (N-1)];
        std::copy(v[0].begin(), v[0].begin() + N * (N-1), m);
        bool parity_swap = false;
        
        if (decomp_plu(m, N-1, N, P, &parity_swap) > 0) {
            // matrix is singular; nullity is > 1. return 0 vector.
            return Vec<T,N>();
        }
        
        for (index_t r = N-2; r >= 0; r--) {
            // back substitute.
            for (index_t c = N-1; c > r; c--) {
                o[r] -= o[c] * _MxElem(m,r,c);
            }
            o[r] /= _MxElem(m,r,r);
        }
        
        return o;
    }
    
    
    template <typename T>
    inline Vec<T,3> orthogonal(const Vec<T,3> v[2]) {
        return v[0] ^ v[1];
    }
    
    
    template <typename T>
    inline Vec<T,2> orthogonal(const Vec<T,2> v[1]) {
        return v[0].leftPerpendicular();
    }
    
    
    /**
     * Compute the null space of a vector basis. 
     * 
     * The computed null bases will not necessarily be orthogonal to each other.
     * Use `orthogonalize()` after computing `nullspace()` if an orthogonal basis
     * is needed.
     * 
     * `bases` and `null_basis` may alias each other.
     *
     * If the any bases are linearly dependent, the resultant `null_basis` will be
     * filled with `N - n` zero vectors.
     * 
     * @param bases Array of `n` linearly independent basis vectors.
     * @param n Number of basis vectors in the array.
     * @param null_basis Array with space to receive `N - n` output bases, whose dot 
     * products with the inputs will all be zero.
     */
    template <typename T, index_t N>
    void nullspace(const Vec<T,N> bases[], index_t n, Vec<T,N> null_basis[]) {
        if (n >= N) return;
        if (N - n == 1) {
            null_basis[0] = orthogonal(bases);
            return;
        }
        
        T       m[N * (N-2)]; // <--  N-2 = max # bases
        index_t P[N - 2];     //     (N-1 case handled above)
              T* n0 = null_basis[0].begin();
        const T* b0 = bases[0].begin();
        
        // copy the bases to new storage
        std::copy(b0, b0 + N * n, m);
        // zero-init the results
        std::fill(n0, n0 + N * (N - n), 0);
        
        // each original basis is a row of:
        detail::MxWrap<T,true>  V = {m, n, N};
        // each new null basis is a column of:
        detail::MxWrap<T,false> X = {n0, N, N - n};
        
        // get V to row-echelon form
        bool parity_swap = false;
        if (decomp_plu(m, n, N, P, &parity_swap) > 0) return;
        
        // foreach null space basis
        for (index_t b = 0; b < N - n; b++) {
            X(b + n, b) = 1;
            // back substitute.
            for (index_t r = n - 1; r >= 0; r--) {
                T& x_i =  X(r, b);
                x_i    = -V(r, n + b);
                for (index_t c = r + 1; c < n; c++) {
                    x_i -= V(r, c) * X(c, b);
                }
                x_i /= V(r, r);
            }
        }
        
        return;
    }


    /**
     * Use the Gram-Schmidt process to orthogonalize a set of basis vectors. 
     * The first basis vector will not change. The remaining vectors may be
     * of arbitrary magnitude, but will be mutually orthogonal to each
     * other and to the first vector.
     * 
     * @param basis Set of `n` bases vectors. 
     * @param n Number of basis vectors between 0 and `N` inclusive.
     */
    template <typename T, index_t N>
    void orthogonalize(Vec<T,N> basis[], index_t n) {
        for (index_t i = 1; i < n; i++) {
            for (index_t j = 0; j < i; j++) {
                // modified gram-schmidt uses the intermediate basis
                // for better numerical stability.
                basis[i] -= basis[i].projectOn(basis[j]);
            }
        }
    }
    
    /**
     * Use the Gram-Schmidt process to orthonormalize a set of basis vectors. 
     * The first basis vector will not change direction. All vectors will
     * be made mutually orthogonal and unit length.
     * 
     * @param basis Set of `n` bases vectors. 
     * @param n Number of basis vectors between 0 and `N` inclusive.
     */
    template <typename T, index_t N>
    void orthonormalize(Vec<T,N> basis[], index_t n) {
        for (index_t i = 1; i < n; i++) {
            for (index_t j = 0; j < i; j++) {
                // modified gram-schmidt uses the intermediate basis
                // for better numerical stability.
                basis[i] -= basis[i].projectOn(basis[j]);
            }
        }
        for (index_t i = 0; i < n; i++) {
            basis[i] /= basis[i].mag();
        }
    }
    
    /// @} //addtogroup linalg
    
} // namespace geom

#undef _MxElem

#endif	/* ORTHOGONAL_H */

