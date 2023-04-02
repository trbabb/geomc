#ifndef _MATRIXDET_H_
#define _MATRIXDET_H_

#include <geomc/Storage.h>
#include <geomc/linalg/LUDecomp.h>

// todo: all untested

namespace geom {
    
namespace detail {

template <typename T>
struct _m2x2 {
    T a,b,
      c,d;
};


template <typename T>
struct _m3x3 {
    T a,b,c,
      d,e,f,
      g,h,i;
};


template <typename T>
struct _m4x4 {
    T a,b,c,d,
      e,f,g,h,
      i,j,k,l,
      m,n,o,p;
};

template <typename T>
inline T cofac(T a, T b, T c, T d, T e, T f) {
    // compute a * b - c * d + e * f
    return geom::multiply_add(e, f, diff_of_products(a, b, c, d));
}

} // namespace detail

/**
 * @addtogroup matrix
 * @{
 */ 

/// Compute the 2x2 matrix determinant from individual parameters.
template <typename T>
inline T det2x2(T a, T b, T c, T d) {
    return diff_of_products(a, d, b, c);
}

template <typename T>
inline T det2x2(const T m[4]) {
    const detail::_m2x2<T>& v = *reinterpret_cast<detail::_m2x2<T>*>(m);
    return diff_of_products(v.a, v.d, v.b, v.c);
}

/// Compute the 3x3 matrix determinant
template <typename T>
inline T det3x3(const T m[9]) {
    const detail::_m3x3<T>& v = *reinterpret_cast<detail::_m3x3<T>*>(m);
    return detail::cofac(
        v.a, det2x2(v.e, v.f, v.h, v.i),
        v.b, det2x2(v.d, v.f, v.g, v.i),
        v.c, det2x2(v.d, v.e, v.g, v.h)
    );
}

/// Compute the 3x3 matrix determinant from individual parameters.
template <typename T>
inline T det3x3(
        T a, T b, T c,
        T d, T e, T f,
        T g, T h, T i)
{
    return detail::cofac(
        a, det2x2(e, f, h, i),
        b, det2x2(d, f, g, i),
        c, det2x2(d, e, g, h)
    );
}

/// Compute the 4x4 matrix determinant.
template <typename T>
T det4x4(const T m[16]) {
    const detail::_m4x4<T>& v = *reinterpret_cast<detail::_m4x4<T>*>(m);
    return 
        diff_of_products(
            v.a, det3x3(v.f, v.g, v.h,  v.j, v.k, v.l,  v.n, v.o, v.p),
            v.b, det3x3(v.e, v.g, v.h,  v.i, v.k, v.l,  v.m, v.o, v.p)
        )
      + diff_of_products(
            v.c, det3x3(v.e, v.f, v.h,  v.i, v.j, v.l,  v.m, v.n, v.p),
            v.d, det3x3(v.e, v.f, v.g,  v.i, v.j, v.k,  v.m, v.n, v.o)
        );
}



/**
 * @brief Destructively compute the determinant of a square matrix.
 * 
 * The provided matrix will be used as a buffer for the computation, and no additional
 * memory will be allocated. The contents of `m` after the function returns are undefined.
 * 
 * Note that it is generally preferable to use one of the `detNxN` methods if the
 * dimension is small and known at compile time.
 * 
 * @param m The square matrix whose determinant is to be computed. `m` may be either row-
 * or column- major layout. The contents of `m` will be overwritten.
 * @param n The number of rows/columns in the matrix.
 */
template <typename T>
T det_destructive(T* m, index_t n) {
    if (n == 1) return m[0];
    bool odd;
    if (decomp_plu(m, n, n, nullptr, &odd) != 0) return 0;
    T d = odd ? -1 : 1;
    // det(A) is the product of its diagonals, if A is triangular
    for (index_t i = 0; i < n; ++i) {
        d *= m[i * i];
    }
    return d;
}

/**
 * @brief Compute the determinant of a square matrix.
 * 
 * If `n` is larger than 8, heap memory may be allocated as buffer space for the computation.
 * (To avoid heap allocation, use `det_destructive()` instead).
 * 
 * @param m The square matrix whose determinant is to be computed. `m` may be either row-
 * or column- major layout.
 * @param n The number of rows/columns in the matrix.
 */
template <typename T>
T detNxN(const T* m, index_t n) {
    bool odd;
    SmallStorage<T, 64> buf(m, n * n);
    if (decomp_plu(buf.get(), n, n, nullptr, &odd) != 0) return 0;
    T d = odd ? -1 : 1;
    // det(A) is the product of its diagonal elems, if A is triangular
    // and det(AB) = det(A) * det(B). So if M = LU,
    // because `L` has unit diagonal, det(M) is just the diagonal elements of `U`.
    for (index_t i = 0; i < n; ++i) {
        d *= buf[i * i];
    }
    return d;
}

/**
 * @brief Compute the determinant of a square matrix.
 * 
 * If the dimension of the matrix is larger than 8, heap memory may be allocated as buffer
 * space for the computation.
 * 
 * If the dimension of the matrix can be determined at compile time to be nonsquare,
 * a static assertion is raised. Otherwise, nonsquare matrices will have determinant 0.
 * 
 * @param m The square matrix whose determinant is to be computed.
 */
template <typename T, index_t M, index_t N, MatrixLayout Lyt, StoragePolicy P>
inline T det(const SimpleMatrix<T, M, N, Lyt, P>& m) {
    static_assert(N * M == 0 or M == N, "Determinant only defined for square matrices.");
    if (M * N == 0 and m.rows() != m.cols()) return 0;
    return det(m.data_begin(), m.rows());
}

// 2x2 matrix det
template <typename T, MatrixLayout Lyt, StoragePolicy P>
inline T det(const SimpleMatrix<T, 2, 2, Lyt, P>& m) {
    const T* v = m.data_begin();
    // n.b.: layout is irrelevant
    return det2x2(v[0], v[1], v[2], v[3]);
}

// 3x3 matrix det
template <typename T, MatrixLayout Lyt, StoragePolicy P>
inline T det(const SimpleMatrix<T, 3, 3, Lyt, P>& m) {
    return det3x3(m.data_begin());
}

// 4x4 matrix det
template <typename T, MatrixLayout Lyt, StoragePolicy P>
inline T det(const SimpleMatrix<T, 4, 4, Lyt, P>& m) {
    return det4x4(m.data_begin());
}

/**
 * @brief Compute the determinant of a diagonal matrix.
 * 
 * If the dimension of the matrix can be determined at compile time to be nonsquare,
 * a static assertion is raised. Otherwise, nonsquare matrices will have determinant 0.
 * 
 * @param m The matrix whose determinant is to be computed.
 */
template <typename T, index_t N, index_t M>
T det(const DiagMatrix<T,N,M>& m) {
    static_assert(N * M == 0 or M == N, "Determinant only defined for square matrices.");
    T d = 1;
    if (m.rows() != m.cols()) return 0;
    for (T* a = m.diagonal_begin(); a != m.diagonal_end(); ++a) {
        d *= *a;
    }
    return d;
}

/**
 * @brief Compute the determinant of a permutation matrix.
 * 
 * If the dimension of the matrix is dynamic and larger than 32, heap memory may
 * be allocated as buffer space for the computation.
 * 
 * If it is desired that no heap memory be allocated, consider using `permutation_sign()`,
 * which is equal to the permutation's determinant, and uses the permutation map
 * (destructively) as a buffer for the calculation.
 * 
 * @param m The permutation matrix whose determinant is to be computed.
 */
template <index_t N>
inline index_t det(const PermutationMatrix<N>& m) {
    index_t n = m.rows();
    SmallStorage<index_t,(N > 32 ? N : 32)> p(m.getRowSources(), n);
    return permutation_sign(p.get(), n);
}

/// @} // ingroup matrix

} // namespace geom


#endif