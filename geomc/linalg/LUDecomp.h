/*
 * LUDecomp.h
 *
 *  Created on: Jul 11, 2013
 *      Author: tbabb
 */

#ifndef LUDECOMP_H_
#define LUDECOMP_H_

#include <type_traits>
#include <geomc/SmallStorage.h>
#include <geomc/linalg/LinalgTypes.h>
#include <geomc/linalg/Vec.h>

#include <geomc/linalg/mtxtypes/SimpleMatrix.h>
#include <geomc/linalg/mtxdetail/MatrixGlue.h>

namespace geom {


// todo: it probably doesn't make sense to allow nonsquare solving of LU-
//       either you want the terms on the other side, or you just get
//       a bunch of zero-rows.
// todo: it probably doesn't make sense to have independent methods for
//       PLU and LUP. If you really want to pivot by columns, then transpose.

namespace detail {

// xxx: todo: can apply non-destructively by using sign bit as a tag
//      but for some reason, it isn't working correctly
// >>> really this should just be implemented/fixed in PermutationMatrix

// apply the permutation map `p` to the n x m matrix `A`, overwriting `p` in the process
// (after which p will be the identity permutation)
template <typename T, bool RowMajor=true>
void destructive_apply_permutation(index_t* p, T* A, index_t n, index_t m) {
    detail::MxWrap<T,RowMajor> mx = {A, n, m};
    for (index_t i = 0; i < n; i++) {
        index_t cur  = i;
        index_t next = p[i]; // source of new element at this position
        if (next == i) continue; // no swap necessary
        // if (next < 0) {
        //     // marked as visited. restore original value
        //     p[i] = -(p[i]);
        //     continue;
        // }
        do {
            for (index_t j = 0; j < m; ++j) {
                // give v[cur] the correct element, v[next] gets the cycle's first element
                std::swap(mx(cur, j), mx(next, j));
            }
            // ...and then proceed to find out what goes into the hole
            cur  = next;
            next = p[next];
            // mark this element as visited by mapping it to itself
            // (we only need to do this if we might encounter the cycle again)
            p[cur] = cur;
            
            // mark as visited
            // if (cur > i) {
            //     p[cur] = -(p[cur]);
            // }
        } while (next != i);
    }
}

} // namespace detail


/** 
 * @ingroup linalg
 * @{
 */


/************************************************
 * Matrix LU decomposition                      *
 *                                              *
 * Implementation using Doolittle's algorithm   *
 ************************************************/

/**
 * @brief LU decomposition, pivoting columns.
 *
 * Solve `LU = MP` for the lower- and upper-triangular matrices L and U
 * and a permutation matrix P. The diagonal of the decomposed matrix
 * belongs to `U` and has arbitrary elements. The diagonals of `L` are
 * implicitly ones.
 *
 * @tparam T The element type of the matrix.
 * @tparam RowMajor Whether the layout of the matrix is row-major 
 * (`true`) or column-major (`false`).
 * 
 * @param m Array of elements to decompose.
 * @param rows Number of rows in `m`.
 * @param cols Number of columns in `m`.
 * @param reorder Array with space for `cols` elements to be filled with the column 
 * source indexes. May be `null`.
 * @param parity Whether an odd number of column-swaps was performed.
 * @return The number of degenerate columns discovered.
 */
template <typename T, bool RowMajor=true>
index_t decomp_lup(T* m, index_t rows, index_t cols, index_t* reorder, bool* parity) {
    const index_t n = std::min(rows, cols);
    index_t degenerate_ct = 0;
    bool _parity = false;
    if (parity == nullptr) parity = &_parity;
    *parity = false;
    // fill permutation array
    if (reorder != nullptr) {
        for (index_t i = 0; i < cols; i++){
            reorder[i] = i;
        }
    }
    
    detail::MxWrap<T,RowMajor> mx = {m, rows, cols};
    
    for (index_t i = 0; i < n - 1; i++) {
        // find pivot
        T biggest = std::abs(mx.elem(i,i));
        index_t pvt = i;
        for (index_t p = i + 1; p < cols; p++) {
            T pvt_val = std::abs(mx.elem(i, p));
            if (pvt_val > biggest) {
                pvt = p;
                biggest = pvt_val;
            }
        }
        
        if (biggest == 0) {
            // singular matrix
            // could test against an epsilon to 
            // find ill-conditioned matrices
            degenerate_ct += 1;
            continue;
        } else if (pvt != i) {
            // swap col <i> with col <pvt>
            for (index_t r = 0; r < rows; r++) {
                std::swap(mx.elem(r, i),
                          mx.elem(r, pvt));
            }
            // make a note of the permutation in <P>
            if (reorder != nullptr) std::swap(reorder[i], reorder[pvt]);
            *parity = !*parity;
        }
        
        // eliminate lower elements
        T a = mx.elem(i,i);
        for (index_t r = i + 1; r < rows; r++) {
            T b = mx.elem(r,i) / a;
            if (b != 0) {
                for (index_t c = i + 1; c < cols; c++) {
                    // R_r = R_r - b * R_i
                    T src_elem = mx.elem(i,c);
                    T dst_elem = mx.elem(r,c);
                    // mx[r,c]   = dst_elem - b * src_elem;
                    mx.elem(r,c) = geom::multiply_add(-b, src_elem, dst_elem);
                }
            }
            // set the lower matrix
            mx.elem(r,i) = b;
        }
    }
    
    return degenerate_ct;
}

/**
 * @brief LU decomposition, pivoting rows.
 * 
 * Solve `LU = PM` for the lower- and upper-triangular matrices L and U
 * and a permutation matrix P. The diagonal of the decomposed matrix
 * belongs to `U` and has arbitrary elements. The diagonals of `L` are
 * implicitly ones.
 *
 * @tparam T The element type of the matrix.
 * @tparam RowMajor Whether the layout of the matrix is row-major 
 * (`true`) or column-major (`false`).
 * 
 * @param m Array of elements to decompose.
 * @param rows Number of rows in `m`.
 * @param cols Number of columns in `m`.
 * @param reorder Array with space for `rows` elements to be filled with
 * the row source indexes. May be `null`.
 * @param parity Whether an odd number of row-swaps was performed.
 * @return The number of degenerate rows discovered.
 */
template <typename T, bool RowMajor=true>
index_t decomp_plu(T* m, index_t rows, index_t cols, index_t* reorder, bool* parity) {
    const index_t n = std::min(rows, cols);
    bool _parity = false;
    if (parity == nullptr) parity = &_parity;
    index_t degenerate_ct = 0;
    *parity = false;
    // fill permutation array
    if (reorder) {
        for (index_t i = 0; i < rows; i++){
            reorder[i] = i;
        }
    }
    
    detail::MxWrap<T,RowMajor> mx = {m, rows, cols};
    
    for (index_t i = 0; i < n - 1; i++) {
        // find pivot
        T biggest = std::abs(mx.elem(i,i));
        index_t pvt = i;
        for (index_t p = i + 1; p < rows; p++) {
            T pvt_val = std::abs(mx.elem(p,i));
            if (pvt_val > biggest) {
                pvt = p;
                biggest = pvt_val;
            }
        }
        
        if (biggest == 0) {
            // singular matrix. could test against an epsilon to 
            // find ill-conditioned matrices
            degenerate_ct += 1;
            continue;
        } else if (pvt != i) {
            // swap row <i> with row <pvt>
            for (index_t c = 0; c < cols; c++) {
                std::swap(mx.elem(i,   c),
                          mx.elem(pvt, c));
            }
            // make a note of the permutation in <P>
            if (reorder) std::swap(reorder[i], reorder[pvt]);
            *parity = !*parity;
        }
        
        // eliminate lower elements
        T a = mx.elem(i,i);
        for (index_t r = i + 1; r < rows; r++) {
            T b = mx.elem(r,i) / a;
            if (b != 0) {
                for (index_t c = i + 1; c < cols; c++) {
                    // R_r = R_r - b * R_i
                    T src_elem = mx.elem(i,c);
                    T dst_elem = mx.elem(r,c);
                    // mx[r,c]   = dst_elem - b * src_elem;
                    mx.elem(r,c) = geom::multiply_add(-b, src_elem, dst_elem);
                }
            }
            // set the lower matrix
            mx.elem(r,i) = b;
        }
    }
    
    return degenerate_ct;
}


/**
 * @brief Solve the decomposed matrix equation `LUX = PB` for `X`, given `LU` and `B`.
 *
 * @tparam T The element type of the matrix.
 * @tparam RowMajor Whether the layout of the matrix is row-major (`true`) 
 * or column-major (`false`).
 *
 * @param plu An `n × n` PLU-decomposed matrix. 
 * @param p The permutation array of row-sources filled by `decomp_plu()`.
 * @param n The number of rows and columns in the matrix.
 * @param m The number of solution columns.
 * @param x The solution matrix of `n × m` elements, having the same row- or
 *   column-major layout as `m`.
 * @param b A matrix of `n × m` elements, having the same layout as `x`.
 * @param skip How many rows in `X`, in order from the first, to skip solving for. 
 * If greater than 0, the corresponding rows within `X` will contain nonsense values.
 */
template <typename T, bool RowMajor=true>
void backsolve_plu(
        const T*       plu,
        const index_t* p,
        index_t        n,
        index_t        m,
        T*             x,
        const T*       b,
        index_t        skip=0)
{
    detail::MxWrap<const T, RowMajor> mx = {plu, n, n};
    detail::MxWrap<const T, RowMajor>  B = {b,   n, m};
    // <y> and <x>'s elements are used such that
    // y[i] is never read after x[i] is written.
    // thus we may collapse their storage and save space:
    detail::MxWrap<T, RowMajor> Y = {x, n, m};
    detail::MxWrap<T, RowMajor> X = {x, n, m};
    
    // LUX = Pb
    // (UX) is a matrix (call it "Y")
    // LY  = Pb
    // so let's solve for each column `i` of Y:
    for (index_t r = 0; r < n; r++) {
        for (index_t i = 0; i < m; ++i) {
            Y.elem(r,i) = B.elem(p[r], i);  // y[r] <- Pb[r]
            for (index_t c = 0; c < r; c++) {
                Y.elem(r,i) -= Y.elem(c,i) * mx.elem(r,c);
            }
        }
    }
    // now with y, we may obtain x from:
    // Ux = y
    for (index_t r = n - 1; r >= skip; r--) {
        // x[r] = y[r]; // nop; x and y are the same!
        T k = mx.elem(r,r);
        for (index_t i = 0; i < m; ++i) {
            for (index_t c = n - 1; c > r; c--) {
                X.elem(r,i) -= X.elem(c,i) * mx.elem(r,c);
            }
            X.elem(r, i) /= k;
        }
    }
}

/**
 * @brief Solve a decomposed system of linear equations `LUX = B`, without a permutation map.
 * 
 * `X` and `B` are `n × m` matrices.
 *
 * @tparam T The element type of the matrix.
 * @tparam RowMajor Whether the layout of the matrix is row-major (`true`) 
 * or column-major (`false`).
 *
 * @param lu An `n × n` LU-decomposed matrix.
 * @param n The number of rows and columns in the matrix.
 * @param m The number of columns in `x` to solve for.
 * @param x The solution vector of `n` elements. This must be initialized to the value of `B`, and
 * will be rewritten to contain the solution `x`.
 * @param skip How many variables, in order from the first, to skip solving for. If greater
 * than 0, the corresponding variables within `x` will contain nonsense values.
 */
template <typename T, bool RowMajor=true>
void backsolve_lu(
        const T* lu,
        index_t  n,
        index_t  m,
        T*       x,
        index_t  skip=0)
{
    detail::MxWrap<const T, RowMajor> mx = {lu, n, n};
    detail::MxWrap<      T, RowMajor> X  = {x,  n, m}; // `x` is initialized with `b`
    
    // LUx = b
    // (Ux) is a vector (call it "y")
    // Ly  = b
    for (index_t r = 1; r < n; r++) {
        // do this for each column i in x:
        for (index_t i = 0; i < m; ++i) {
            // x[r,i] = b[r,i]  // <-- already true
            for (index_t c = 0; c < r; c++) {
                X.elem(r,i) -= X.elem(c,i) * mx.elem(r,c);
            }
        }
    }
    // now with y, we may obtain x from:
    // Ux = y
    for (index_t r = n - 1; r >= skip; r--) {
        T k = mx.elem(r,r);
        // do for each column i in y
        for (index_t i = 0; i < m; ++i) {
            for (index_t c = n - 1; c > r; c--) {
                X.elem(r,i) -= X.elem(c,i) * mx.elem(r,c);
            }
            X.elem(r, i) /= k;
        }
    }
}


/**
 * @brief Solve the matrix equation `AX = B` for the matrix `X`, given matrices `A` and `B`.
 *
 * @tparam T The element type of the matrix.
 * @tparam RowMajor Whether the layout of the matrix is row-major (`true`)
 * or column-major (`false`).
 * 
 * @param a A buffer of elements in the n × n matrix `A`. This array will
 * be altered during the solution process, so pass a copy if the original 
 * needs to remain unchanged.
 * @param n The number of rows in the matrix.
 * @param m The number of columns in `X` and `B`.
 * @param x The matrix of `n × m` initially containing `b`, which is to be filled with
 * the solution vector `x`.
 * @param skip How many rows in `X`, in order from the first, to skip solving for. If greater 
 * than 0, the corresponding rows within `X` will contain nonsense values.
 * @return `true` if a solution could be found; `false` if the system is degenerate.
 */
template <typename T, bool RowMajor=true>
inline bool linear_solve(T* a, index_t n, index_t m, T* x, index_t skip=0) {
    SmallStorage<index_t, 24> p(n); // unlikely to need an alloc
    bool parity;
    if (decomp_plu<T,RowMajor>(a, n, n, p.begin(), &parity) > 0) {
        return false;
    }
    // pre-apply P to x
    detail::destructive_apply_permutation<T,RowMajor>(p.begin(), x, n, m);
    // LUx = Pb
    backsolve_lu<T,RowMajor>(a, n, m, x, skip);
    return true;
}


/**
 * @brief Write each vector in `b` in terms of the basis vectors in `bases`.
 * 
 * Return an `x_j` for each such `b_j` that `sum(bases[i] * x_j[i]) = b_j`.
 * 
 * @param bases An array of `N` basis vectors. The contents of this array will
 * be altered during the solution process, so pass a copy if the original 
 * array needs to remain unchanged.
 * @param m The number of vectors in `x` to solve for.
 * @param x An array of `m` vectors `b_j`, to be overwritten with the `x_j` solutions.
 * @param skip How many rows, in order from the first, to skip solving for. If greater
 * than 0, the corresponding variables within each `x` will contain nonsense values.
 * @return `true` if a solution could be found; `false` if the system is degenerate.
 */
template <typename T, index_t N>
inline bool linear_solve(Vec<T,N> bases[N], index_t m, Vec<T,N>* x, index_t skip=0) {
    T* _x = x->begin();
    if (N < 5) {
        // matrix inv is empirically faster than solve() for N < 5.
        WrapperMatrix<T,N,N,COL_MAJOR> mx(bases->begin());
        WrapperMatrix<T,N,0,COL_MAJOR>  X(_x, m);
        SimpleMatrix<T,N,N> m_inv;
        if (not inv(&m_inv, mx)) return false;
        mul(&X, m_inv, X);
        return true;
    } else {
        index_t p[N];
        T* const mx = bases[0].begin();
        bool _parity;
        // LU = PM (M's rows reordered)
        if (decomp_plu<T,false>(mx, N, N, p, &_parity) > 0) {
            return false;
        }
        // pre-apply P to x
        detail::destructive_apply_permutation<T,false>(p, _x, N, m);
        // LUx = Pb
        backsolve_lu<T,false>(mx, N, m, _x, skip);
        return true;
    }
}


/// @}  ingroup linalg

}  // namespace geom

#endif /* LUDECOMP_H_ */
