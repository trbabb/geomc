#pragma once
/*
 * MatrixInv.h
 *
 *  Created on: Jul 18, 2013
 *      Author: tbabb
 */

#include <algorithm>

#include <geomc/SmallStorage.h>
#include <geomc/linalg/LinalgTypes.h>
#include <geomc/linalg/LUDecomp.h>
#include <geomc/linalg/mtxdetail/MatrixDet.h>
#include <geomc/function/Utils.h>

namespace geom {


/*************************************************
 * General matrix inversion                      *
 *                                               *
 * Inverts the data in an iterator-like object   *
 *************************************************/

// todo: for N=3,4 and T in {float,double}, we could implement SSE-optimized 
//       versions of the matrix inverters. Unclear if that would beat -O3.
// todo: all of this is untested.
// todo: it's not currently possible to amortize the allocation of P
//       for large matrices.

/**
 * @addtogroup matrix
 * @{
 */

/**
 * @brief 2 × 2 matrix inversion on a flat array.
 * 
 * `out` and `m` may alias each other. It is assumed that `m` and `out` have the same
 * layout.
 * 
 * @param out An array with space for four elements to receive the computed inverse.
 * @param m The 2 × 2 matrix to invert.
 * @return `true` if the matrix is invertible, false otherwise. If the matrix is not
 * invertible, then `out` is unchanged.
 */
template <typename T>
bool inv2x2(T* out, const T m[4]) {
    T a = m[0];
    T b = m[1];
    T c = m[2];
    T d = m[3];
    T det = det2x2(a, b, c, d);
    if (det == 0) return false;
    T inv[4] = { d / det, -b / det,
                -c / det,  a / det };
    std::copy(inv, inv + 4, out);
    return true;
}


/**
 * @brief 3 × 3 matrix inversion on a flat array.
 * 
 * `out` and `m` may alias each other. It is assumed that `m` and `out` have the
 * same layout.
 * 
 * @param out An array with space for nine elements to receive the computed inverse.
 * @param m The 3 × 3 matrix to invert.
 * @return `true` if the matrix is invertible, false otherwise. If the matrix is not
 * invertible, then `out` is unchanged.
 */
template <typename T>
bool inv3x3(T* out, const T m[9]) {
    const detail::_m3x3<T>& s = *reinterpret_cast<const detail::_m3x3<T>*>(m);
    // inverse 3x3 matrix by cramer's rule
    T p0 = det2x2(s.i, s.h, s.f, s.e);
    T p1 = det2x2(s.i, s.h, s.c, s.b);
    T p2 = det2x2(s.f, s.e, s.c, s.b);
    
    // det = det2x2(s.a, s.d, p1, p0) + s.g * p2;
    T det = geom::multiply_add(s.g, p2, det2x2(s.a, s.d, p1, p0));
    if (det == 0) return false;
    
    T inv[9] = {
        det2x2(s.i, s.h, s.f, s.e) / det,
        det2x2(s.h, s.i, s.b, s.c) / det,
        det2x2(s.f, s.e, s.c, s.b) / det,
    
        det2x2(s.g, s.i, s.d, s.f) / det,
        det2x2(s.i, s.g, s.c, s.a) / det,
        det2x2(s.d, s.f, s.a, s.c) / det,
    
        det2x2(s.h, s.g, s.e, s.d) / det,
        det2x2(s.g, s.h, s.a, s.b) / det,
        det2x2(s.e, s.d, s.b, s.a) / det
    };
    
    std::copy(inv, inv + 9, out);
    return true;
}


/**
 * @brief 4 × 4 matrix inversion on a flat array.
 * 
 * `out` and `m` may alias each other. It is assumed that `m` and `out` have the
 * same layout.
 * 
 * @param out An array with space for 16 elements to receive the computed inverse.
 * @param m The 4 × 4 matrix to invert.
 * @return `true` if the matrix is invertible, false otherwise. If the matrix is not
 * invertible, then `out` is unchanged.
 */
template <typename T>
bool inv4x4(T* out, const T m[16]) {
    const detail::_m4x4<T>& s = *reinterpret_cast<const detail::_m4x4<T>*>(m);
    T s0 = det2x2(s.a, s.b, s.e, s.f);
    T s1 = det2x2(s.a, s.c, s.e, s.g);
    T s2 = det2x2(s.a, s.d, s.e, s.h);
    T s3 = det2x2(s.b, s.c, s.f, s.g);
    T s4 = det2x2(s.b, s.d, s.f, s.h);
    T s5 = det2x2(s.c, s.d, s.g, s.h);
    
    T c0 = det2x2(s.i, s.j, s.m, s.n);
    T c1 = det2x2(s.i, s.k, s.m, s.o);
    T c2 = det2x2(s.i, s.l, s.m, s.p);
    T c3 = det2x2(s.j, s.k, s.n, s.o);
    T c4 = det2x2(s.j, s.l, s.n, s.p);
    T c5 = det2x2(s.k, s.l, s.o, s.p);
    
    // det = s0*c5 - s1*c4 + s2*c3 + s3*c2 - s4*c1 + s5*c0;
    T det = diff_of_products(s0, c5, s1, c4) +
            diff_of_products(s3, c2, s4, c1) +
             sum_of_products(s2, c3, s5, c0);
    
    if (det == 0) return false;
    
    T inv[16] = {
            // cofac is:   a * b  -  c * d  +  e * f
            detail::cofac(s.f, c5,  s.g, c4,  s.h, c3) / det,
           -detail::cofac(s.b, c5,  s.c, c4,  s.d, c3) / det,
            detail::cofac(s.n, s5,  s.o, s4,  s.p, s3) / det,
           -detail::cofac(s.j, s5,  s.k, s4,  s.l, s3) / det,
             
           -detail::cofac(s.e, c5,  s.g, c2,  s.h, c1) / det,
            detail::cofac(s.a, c5,  s.c, c2,  s.d, c1) / det,
           -detail::cofac(s.m, s5,  s.o, s2,  s.p, s1) / det,
            detail::cofac(s.i, s5,  s.k, s2,  s.l, s1) / det,
             
            detail::cofac(s.e, c4,  s.f, c2,  s.h, c0) / det,
           -detail::cofac(s.a, c4,  s.b, c2,  s.d, c0) / det,
            detail::cofac(s.m, s4,  s.n, s2,  s.p, s0) / det,
           -detail::cofac(s.i, s4,  s.j, s2,  s.l, s0) / det,
             
           -detail::cofac(s.e, c3,  s.f, c1,  s.g, c0) / det,
            detail::cofac(s.a, c3,  s.b, c1,  s.c, c0) / det,
           -detail::cofac(s.m, s3,  s.n, s1,  s.o, s0) / det,
            detail::cofac(s.i, s3,  s.j, s1,  s.k, s0) / det
    };
    
    std::copy(inv, inv + 16, out);
    return true;
}


/**
 * @brief N × N matrix inversion on a flat array.
 * 
 * After the inversion, `m` will have undefined contents, regardless of
 * whether it was invertible.
 * 
 * It is assumed that `m` and `out` have the same layout. `m` *must not* alias `out`.
 * 
 * This function uses a general inversion algorithm, and the specialized 
 * dimension-specific algorithms will be preferable when the dimension is known.
 * 
 * @param out An array with space for `N`<sup>2</sup> elements to receive the computed inverse.
 * @param m The N × N matrix to invert.
 * @param n The dimension `N` of the matrix.
 * @return `true` if the matrix is invertible, false otherwise. If the matrix is not
 * invertible, then `out` is unchanged.
 */
template <typename T>
bool invNxN(T* out, T* m, index_t n) {
    SmallStorage<index_t, 24> p(n);
    bool parity = false;
    if (decomp_plu(m, n, n, p.begin(), &parity) > 0) return false;
    std::fill(out, out + n * n, (T)0);
    // pre-apply the permutation P to the identity matrix:
    detail::MxWrap<T,true> o = {out, n, n};
    for (index_t i = 0; i < n; ++i) {
        o(i, p[i]) = 1;
    }
    backsolve_lu(m, n, n, out);
    return true;
}


/**
 * @brief N × N matrix inversion.
 * 
 * After the inversion, `m` will have undefined contents, regardless of
 * whether it was invertible. 
 * 
 * It is assumed that `m` and `out` have the same layout. `m` *must not* alias `out`.
 * 
 * This function uses a general inversion algorithm, and the specialized 
 * dimension-specific algorithms will be preferable when the dimension is known.
 * 
 * @tparam N The dimension `N` of the matrix.
 * @param out An array with space for `N`<sup>2</sup> elements to receive the computed inverse.
 * @param m The N × N matrix to invert.
 * @return `true` if the matrix is invertible, false otherwise. If the matrix is not
 * invertible, then `out` is unchanged.
 */
template <typename T, index_t N>
bool invNxN(T* out, T* m) {
    index_t p[N];
    bool parity = false;
    if (decomp_plu(m, N, N, p, &parity) > 0) return false;
    std::fill(out, out + N * N, (T)0);
    // pre-apply the permutation P to the identity matrix:
    detail::MxWrap<T,true> o = {out, N, N};
    for (index_t i = 0; i < N; ++i) {
        o(i, p[i]) = 1;
    }
    backsolve_lu<T,true>(m, N, N, out);
    return true;
}

/// @} // addtogroup matrix


namespace detail {


// NxN static inverse spec
template <typename T, index_t N>
struct _ImplMtxInvRaw {
    template <typename Md, typename Mx>
    static inline bool inv(Md* out, const Mx& m) {
        // copy to a buffer having the same layout as the destination matrix:
        SimpleMatrix<T,N,N,Md::Layout> buf(out->rows(), out->rows());
        detail::_mtxcopy(&buf, m);
        // do the inverse, destroying the buffer
        return geom::invNxN<T,N>(out->data_begin(), buf.data_begin());
    }
};


// dynamic matrix inverse
template <typename T>
struct _ImplMtxInvRaw<T,0> {
    
    template <typename Md, typename Mx>
    static bool inv(Md* into, const Mx& m) {
        const index_t n = m.rows();
        switch(n) {
            case 2:
                return geom::detail::_ImplMtxInvRaw<T,2>::inv(into, m);
            case 3:
                return geom::detail::_ImplMtxInvRaw<T,3>::inv(into, m);
            case 4:
                return geom::detail::_ImplMtxInvRaw<T,4>::inv(into, m);
            default:
                // make a buffer with the same layout as the output
                static constexpr MatrixLayout Lyt = Md::Layout;
                SmallStorage<T,64> buf(n * n); // try not to alloc
                WrapperMatrix<T,0,0,Lyt> bmx(buf.get(), n, n);
                detail::_mtxcopy(&bmx, m);
                return geom::invNxN(into->data_begin(), bmx.data_begin(), n);
        }
    }
};


// 4x4 inverse spec
template <typename T>
struct _ImplMtxInvRaw<T,4> {
    
    template <
        index_t M1, index_t N1, MatrixLayout Lyt1, StoragePolicy P1,
        index_t M2, index_t N2, MatrixLayout Lyt2, StoragePolicy P2>
    static bool inv(
                  SimpleMatrix<T,M1,N1,Lyt1,P1>* out, 
            const SimpleMatrix<T,M2,N2,Lyt2,P2>& m)
    {
        if (not geom::inv4x4(out->data_begin(), m.data_begin())) return false;
        if (Lyt1 != Lyt2) {
            geom::transpose_square_matrix<T,4>(out->data_begin());
        }
        return true;
    }
    
};


// 3x3 inverse spec
template <typename T>
struct _ImplMtxInvRaw<T,3> {
    
    template <
        index_t M1, index_t N1, MatrixLayout Lyt1, StoragePolicy P1,
        index_t M2, index_t N2, MatrixLayout Lyt2, StoragePolicy P2>
    static bool inv(
                  SimpleMatrix<T,M1,N1,Lyt1,P1>* out,
            const SimpleMatrix<T,M2,N2,Lyt2,P2>& m)
    {
        if (not geom::inv3x3(out->data_begin(), m.data_begin())) return false;
        if (Lyt1 != Lyt2) {
            geom::transpose_square_matrix<T,3>(out->data_begin());
        }
        return true;
    }
};


// 2x2 inverse spec
template <typename T>
struct _ImplMtxInvRaw<T,2>{
    template <
        index_t M1, index_t N1, MatrixLayout Lyt1, StoragePolicy P1,
        index_t M2, index_t N2, MatrixLayout Lyt2, StoragePolicy P2>
    static bool inv(
                  SimpleMatrix<T,M1,N1,Lyt1,P1>* out,
            const SimpleMatrix<T,M2,N2,Lyt2,P2>& m)
    {
        if (not geom::inv2x2(out->data_begin(), m.data_begin())) return false;
        if (Lyt1 != Lyt2) {
            T* v = out->data_begin();
            std::swap(v[1], v[2]);
        }
        return true;
    }
};


/*************************************
 * Matrix class inversion            *
 *************************************/

//TODO: for icc compiler with floats, use:
//        ftp://download.intel.com/design/PentiumIII/sml/24504301.pdf
//      or:
//        https://github.com/LiraNuna/glsl-sse2/blob/master/source/mat4.h#L324

template <typename Mx>
class _ImplMtxInv {
public:
    
    static constexpr index_t N = std::max(Mx::ROWDIM, Mx::COLDIM);
    typedef typename Mx::elem_t T;
    typedef geom::SimpleMatrix<T,N,N> return_t;
    
    template <typename Md>
    static inline bool inv(Md* into, const Mx& m) {
        // if any matrix has a static dimension, assume it's been verified to be
        // the correct dimension for the operation, and use that as the static
        // dimension for everybody:
        static constexpr index_t K = std::max(
            std::max(Md::ROWDIM, Md::COLDIM),
            N);
        return geom::detail::_ImplMtxInvRaw<T,K>::inv(into, m);
    }
};

// diagonal matrices are trivial to invert
template <typename T, index_t M, index_t N>
class _ImplMtxInv<geom::DiagMatrix<T,M,N> > {
public:
    
    typedef geom::DiagMatrix<T, (M != 0 ? M : N), (M != 0 ? M : N)> return_t;
    
    template <typename Md>
    static bool inv(Md *into, const geom::DiagMatrix<T,M,N>& m) {
        into->set_zero(); // for newly allocated matrices, redundant. oh well.
        T *p = m.diagonal_begin();
        for (index_t i = 0; p != m.diagonal_end(); ++i, ++p) {
            T elem = *p;
            if (elem == 0) return false;
            into->set(i,i, 1/elem);
        }
        return true;
    }
};

// inverse of permutation matrix is its transpose.
template <index_t N>
class _ImplMtxInv<geom::PermutationMatrix<N>> {
    typedef geom::PermutationMatrix<N> M;
public:
    
    typedef geom::PermutationMatrix<N> return_t;
    
    template <typename Md>
    static inline bool inv(Md *into, const M& m) {
        _ImplMtxTxpose<M>::transpose(into, m);
        return true;
    }
};

// workaround likely clang bug.
// clang tries to make a bogus template substitution (Matrix for bool)
// and chokes.
template <typename T>
class _ImplMtxInv<T*> {};

} // namespace detail
} // namespace geom
