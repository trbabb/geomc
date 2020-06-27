#ifndef __CHOLESKY_H__
#define __CHOLESKY_H__

#include <geomc/linalg/Matrix.h>

// todo: make a Symmetric matrix class, which this can operate on.
//   also optimize addition / multiplication / inversion / det / etc.
//   also author a square(&symm, mtx) function to compute M * M^T into a symm. mtx
//     (and work all of the above into Kalman filter).
// todo: consider also a Triangular matrix class.
// todo: there is no need to template MatrixLayout— square symmetric matrices
//       can be safely transposed.
// todo: use data_begin() and data_end() for the matrix.

namespace geom {

/** @ingroup linalg
 *  @{
 */

/**
 * Perform an in-place Cholesky decomposition on the given square positive-definite
 * matrix `A`, producing a lower-triangular matrix `M` such that 
 * `(M * M`<sup>T</sup>`) = A`.
 *
 * @tparam T Element type.
 * @tparam RowMajor Whether the elements in `m` are arranged in row-major (true) 
 * or column-major (false) order.
 * @param m The elements of the matrix to be decomposed.
 * @param n The number of rows/columns in the matrix to be decomposed.
 * @return True if the decomposition could be completed; false if the matrix was not 
 * positive-definite or ill-conditioned.
 */
template <typename T, bool RowMajor>
bool cholesky(T* m, index_t n) {
    bool ok = true;
    detail::MxWrap<T,RowMajor> mx = {m, n, n};
    for (index_t r = 0; r < n; ++r) {
        for (index_t c = 0; c <= r; ++c) {
            T s = mx.elem(r, c);
            for (index_t k = 0; k < c; ++k) {
                s -= mx.elem(r, k) * mx.elem(c, k);
            }
            if (r == c) {
                ok = ok and (s > 0);
                s = std::sqrt(s);
            } else {
                s = s / mx.elem(c, c);
                mx.elem(c, r) = 0;
            }
            mx.elem(r, c) = s;
        }
    }
    return ok;
}


/**
 * Perform an in-place Cholesky decomposition on the given square positive-definite 
 * matrix `A`, producing a lower-triangular matrix `M` such that 
 * `(M * M`<sup>T</sup>`) = A`.
 * 
 * @param m A square positive-definite matrix.
 * @return False if the matrix is not square, not positive-definite, or could
 * not be decomposed due to ill-conditioning; `true` if the decomposition was
 * completed successfully.
 */
#ifdef PARSING_DOXYGEN
template <typename T, index_t M, index_t N>
bool cholesky(SimpleMatrix<T,M,N>* m);
#else
template <typename T, index_t M, index_t N, MatrixLayout Lyt, StoragePolicy P>
inline typename boost::enable_if_c<M == N or M * N == 0, bool>::type
cholesky(SimpleMatrix<T,M,N,Lyt,P>* m) {
    if (M * N == 0 and m->rows() != m->cols()) {
        return false;
    }
    return cholesky<T,true>(m->begin(), m->rows());
}
#endif


/// @}  // ingroup linalg
    
}

#endif
