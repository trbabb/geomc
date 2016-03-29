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
           1   0

dx + ey + fz + g = 0
0x + hy + iz + j = 0
0x + 0y + kz + l = 0
--
          kz = -l
           z = -l / k
           ...
           etc.
 */

#define _MxElem(m,r,c) m[N*r + c]

namespace geom {
    
    // todo: orthogonal() does not respect winding. the wedge product is maybe something you need.
    
    /**
     * @addtogroup linalg
     * @{
     */
    
    /**
     * Return a vector orthogonal to the given `N-1` vectors.
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
        
        if (decompLUP(m, N-1, N, P, &parity_swap)) {
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
        
        PermutationMatrix<N> Pmtx;
        Pmtx.setRowSources(P);
        
        // invert the permutation.
        // TODO: is that parity_swap correct?
        return (parity_swap ? -1 : 1) * o * Pmtx;
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
     * `bases` and `null_basis` may alias each other.
     * 
     * @param bases Array of `n` linearly independent basis vectors.
     * @param n Number of basis vectors in the array.
     * @param null_basis Array with space for `N - n` output bases, whose dot 
     * products with the inputs are all zero. The elements of this array will not
     * necessarily be orthogonal to each other.
     */
    template <typename T, index_t N>
    void nullspace(const Vec<T,N> bases[], index_t n, Vec<T,N> null_basis[])
    {
        if (n >= N) return;
        if (N - n == 1) {
            null_basis[0] = orthogonal(bases);
            return;
        }
        
        const T *b0 = bases[0].begin();
        index_t P[N];
        T       m0[N*N];
        T * const m1 = m0 + (n*N); // space for new basis
        std::copy(b0, b0 + (n*N),     m0);
        std::fill(m1, m1 + (N-n)*N, (T)0);
        bool parity_swap = false;
        
        decompLUP(m0, n, N, P, &parity_swap);
        
        // todo: handle extra degenerate bases.
        // todo: refactor this to use a mtx wrapper.
        
        // foreach null space basis
        for (index_t b = 0; b < N - n; b++) {
            _MxElem(m1,b,N-b-1) = 1;
            // back substitute.
            for (index_t r = n - 1; r >= 0; r--) {
                for (index_t c_ct = 0; c_ct < n - r; c_ct++) {
                    index_t c = c_ct == 0 ? N - b - 1 : n - c_ct;
                    _MxElem(m1,b,r) -= _MxElem(m1,b,c) * _MxElem(m0,r,c);
                }
                _MxElem(m1,b,r) /= _MxElem(m0,r,r);
            }
        }
        
        // just use this guy to invert the permutation:
        PermutationMatrix<N> Pmtx;
        Pmtx.setRowSources(P);
        const index_t *P_inv = Pmtx.getColSources();
        
        // todo: parity swap ???
        for (index_t b = 0; b < N - n; b++) {
            for (index_t i = 0; i < N; i++) {
                null_basis[b][i] = _MxElem(m1, b, P_inv[i]);
            }
        }
    }
    
    /**
     * Use the Gram-Schmidt process to orthonormalize a set of basis vectors. 
     * The first basis vector will not change direction. 
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

