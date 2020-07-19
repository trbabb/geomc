/*
 * LUDecomp.h
 *
 *  Created on: Jul 11, 2013
 *      Author: tbabb
 */

#ifndef LUDECOMP_H_
#define LUDECOMP_H_

#include <geomc/linalg/Matrix.h>

namespace geom {

// todo: we compute the inverse of P as part of matrix inversion.
//       I suspect in some circumstances (i.e. certain matrix layout
//       combinations) the inverse isn't needed. Can we avoid computing
//       P's inverse?
//       Generally, we should always be able to pre- or post-permute a 
//       vector cheaply, since they are inverses of each other.

// todo: we do a final transpose as part of matrix inversion when
//       the destination is row-major. is there a way we can eliminate
//       this extra op?
//       > note: You can't just transpose LU; that changes its meaning
//         (because multiplication/decomposition does not commute with a transpose)
//       > yes, we might be able to do this, if we break out inverse().
//         we can pre-transpose M before decomposing, by flipping M's layout just
//         before we decompose it. a layout flip can be zero-cost because it
//         is simply a template parameter on the decompose & solve funtions.
//         we'll have both the src and dst matrix in hand, so we'll have enough
//         information do choose the flips correctly. That way, LU will be the 
//         decomposition of the transpose (instead of the erroneous transpose of the 
//         decompostion).
//         I think we must also flip PI (so that we are solving for the ultimate rows 
//         of PI instead of its cols)— we can get this with our choice of P's 
//         internal indexing array.

/** @ingroup linalg
 *  @{
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
 * and a permutation matrix P.
 *
 * @tparam T The element type of the matrix.
 * @tparam RowMajor Whether the layout of the matrix is row-major 
 * (`true`) or column-major (`false`).
 * 
 * @param m Array of elements to decompose.
 * @param rows Number of rows in `m`.
 * @param cols Number of columns in `m`.
 * @param reorder Array with space for `cols` elements to be filled with the column 
 * source indexes.
 * @param swap_parity Whether an odd number of column-swaps was performed.
 * @return The number of degenerate columns discovered.
 */
template <typename T, bool RowMajor=true>
index_t decomp_lup(T* m, index_t rows, index_t cols, index_t* reorder, bool* swap_parity) {
    const index_t n = std::min(rows, cols);
    index_t degenerate_ct = 0;
    *swap_parity = false;
    // fill permutation array
    for (index_t i = 0; i < cols; i++){
        reorder[i] = i;
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
            std::swap(reorder[i], reorder[pvt]);
            *swap_parity = !*swap_parity;
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
                    mx.elem(r,c) = dst_elem - b * src_elem;
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
 * and a permutation matrix P.
 *
 * @tparam T The element type of the matrix.
 * @tparam RowMajor Whether the layout of the matrix is row-major 
 * (`true`) or column-major (`false`).
 * 
 * @param m Array of elements to decompose.
 * @param rows Number of rows in `m`.
 * @param cols Number of columns in `m`.
 * @param reorder Array with space for `rows` elements to be filled with
 * the row source indexes.
 * @param swap_parity Whether an odd number of row-swaps was performed.
 * @return The number of degenerate rows discovered.
 */
template <typename T, bool RowMajor=true>
index_t decomp_plu(T* m, index_t rows, index_t cols, index_t* reorder, bool* swap_parity) {
    const index_t n = std::min(rows, cols);
    index_t degenerate_ct = 0;
    *swap_parity = false;
    // fill permutation array
    for (index_t i = 0; i < rows; i++){
        reorder[i] = i;
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
            std::swap(reorder[i], reorder[pvt]);
            *swap_parity = !*swap_parity;
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
                    mx.elem(r,c) = dst_elem - b * src_elem;
                }
            }
            // set the lower matrix
            mx.elem(r,i) = b;
        }
    }
    
    return degenerate_ct;
}


/**
 * Solve a system of linear equations `LUx = Pb`.
 *
 * @tparam T The element type of the matrix.
 * @tparam RowMajor Whether the layout of the matrix is row-major (`true`) 
 * or column-major (`false`).
 *
 * @param plu An `n x n` PLU-decomposed matrix. 
 * @param p The permutation array of row-sources filled by `decomp_plu()`.
 * @param n The number of rows and columns in the matrix.
 * @param x The solution vector of `n` elements to be filled.
 * @param b A vector of `n` elements.
 * @param skip How many variables, in order from the first, to skip solving for. 
 * If greater than 0, the corresponding variables within `x` will contain nonsense values.
 */
template <typename T, bool RowMajor=true>
void linear_solve_plu(
        const T* plu,
        const index_t* p,
        index_t n,
        T* x,
        const T* b,
        index_t skip=0)
{
    detail::MxWrap<const T, RowMajor> mx = {plu, n, n};
    
    //       LU   = PM
    //       LU x = Pb  <-- our equation.
    // (P^-1)LU x = b
    //        M x = b   <-- same `x`.
    
    // <y> and <x>'s elements are used such that
    // y[i] is never read after x[i] is written.
    // thus we may collapse their storage and save space:
    T* y = x;
    
    // LUx = Pb
    // (Ux) is a vector (call it "y"), so let's solve for it:
    // Ly  = Pb
    y[0] = b[p[0]];
    for (index_t r = 1; r < n; r++) {
        y[r] = b[p[r]];  // y[r] <- Pb[r]
        for (index_t c = 0; c < r; c++) {
            y[r] -= y[c] * mx.elem(r,c);
        }
    }
    // now with y, we may obtain x from:
    // Ux = y
    for (index_t r = n - 1; r >= skip; r--) {
        // x[r] = y[r]; // nop; x and y are the same!
        for (index_t c = n - 1; c > r; c--) {
            x[r] -= x[c] * mx.elem(r,c);
        }
        x[r] /= mx.elem(r,r);
    }
}


/**
 * Solve a system of linear equations `LUx = b`, without a permutation map.
 *
 * @tparam T The element type of the matrix.
 * @tparam RowMajor Whether the layout of the matrix is row-major (`true`) 
 * or column-major (`false`).
 *
 * @param lu An `n x n` LU-decomposed matrix.
 * @param n The number of rows and columns in the matrix.
 * @param x The solution vector of `n` elements to be filled.
 * @param b A vector of `n` elements.
 * @param skip How many variables, in order from the first, to skip solving for. If greater
 * than 0, the corresponding variables within `x` will contain nonsense values.
 */
template <typename T, bool RowMajor=true>
void linear_solve_lu(const T* lu, index_t n, T* x, const T* b, index_t skip=0) {
    detail::MxWrap<const T,RowMajor> mx = {lu, n, n};
    
    // <y> and <x>'s elements are used such that
    // y[i] is never read after x[i] is written.
    // thus we may collapse their storage and save space:
    T* y = x;
    
    // LUx = b
    // (Ux) is a vector, so let's solve for it:
    // Ly  = b
    y[0] = b[0];
    for (index_t r = 1; r < n; r++) {
        y[r] = b[r];
        for (index_t c = 0; c < r; c++) {
            y[r] -= y[c] * mx.elem(r,c);
        }
    }
    // now with y, we may obtain x from:
    // Ux = y
    for (index_t r = n - 1; r >= skip; r--) {
        // x[r] = y[r]; // nop; x and y are the same!
        for (index_t c = n - 1; c > r; c--) {
            x[r] -= x[c] * mx.elem(r,c);
        }
        x[r] /= mx.elem(r,r);
    }
}


/**
 * Solve a system of linear equations `Mx = b`.
 *
 * @tparam T The element type of the matrix.
 * @tparam RowMajor Whether the layout of the matrix is row-major (`true`)
 * or column-major (`false`).
 * 
 * @param m A buffer of elements in the matrix `M`. This array will
 * be altered during the solution process, so pass a copy if the original 
 * needs to remain unchanged.
 * @param n The number of rows in the matrix.
 * @param x The solution vector of `N` elements to be filled.
 * @param b A vector of `N` elements. 
 * @param skip How many variables, in order from the first, to skip solving for. If greater 
 * than 0, the corresponding variables within `x` will contain nonsense values.
 * @return `true` if a solution could be found; `false` if the system is degenerate.
 */
template <typename T, bool RowMajor=true>
inline bool linear_solve(T* m, index_t n, T* x, const T* b, index_t skip=0) {
    SmallStorage<index_t, 24> p(n); // unlikely to need an alloc
    bool parity;
    if (decomp_plu<T,RowMajor>(m, n, n, p.get(), &parity) > 0) {
        return false;
    }
    linear_solve_plu<T,RowMajor>(m, p.get(), n, x, b, skip);
    return true;
}


/**
 * Write a vector `b` in terms of the basis vectors in `bases`; return `x` such that 
 * `sum(bases[i] * x[i]) = b`.
 * 
 * @param bases An array of `N` basis vectors. The contents of this array will
 * be altered during the solution process, so pass a copy if the original 
 * array needs to remain unchanged.
 * @param x The solution vector to be filled.
 * @param b A vector.
 * @param skip How many variables, in order from the first, to skip solving for. If greater
 * than 0, the corresponding variables within `x` will contain nonsense values.
 * @return `true` if a solution could be found; `false` if the system is degenerate.
 */
template <typename T, index_t N>
inline bool linear_solve(Vec<T,N>* bases, Vec<T,N>* x, const Vec<T,N>& b, index_t skip=0) {
    if (N < 5) {
        // matrix inv is empirically faster than solve() for N < 5.
        
        // m is row major, but our bases would be columns.
        // m and its inverse are thus transposed.
        WrapperMatrix<T,N,N> m(bases[0].begin());
        SimpleMatrix<T,N,N>  m_inv_txpose;
        if (!inv(&m_inv_txpose, m)) return false;
        // reverse mult order to get mul by txpose:
        *x = b * m_inv_txpose;
        return true;
    } else {
        index_t p[N];
        T* const m = bases[0].begin();
        bool _parity;
        // LU = PM (M's rows reordered)
        if (decomp_plu<T,false>(m, N, N, p, &_parity) > 0) {
            return false;
        }
        // LUx = Pb
        linear_solve_plu<T,false>(m, p, N, x->begin(), b.begin(), skip);
        
        return true;
    }
}


/// @}  ingroup linalg


//////////// PLU class ////////////

/** 
 * @ingroup matrix
 * @brief Computes the PLU decompostion for a matrix `M`, such that `PM = LU`.
 * 
 * `L` and `U` are lower and upper triangular matrices, respectively.
 * `P` has dimension `(LU.rows() x LU.rows())`, and `LU` has the dimension of `M`.
 */
template <typename T, index_t M, index_t N>
class PLUDecomposition {
protected:
    
    // Dimension of the diagonal of `LU`.
    static const index_t DIAG = (M<N)?M:N;
    // Matrix type for holding `L`.
    typedef SimpleMatrix<T,M,DIAG> L_t;
    // Matrix type for holding `U`.
    typedef SimpleMatrix<T,DIAG,N> U_t;
    
public:

    /**
      * @brief The upper and lower triangular matrices L and U, stored superimposed
      * into each other using a single matrix as storage.
      *
      * The elements of `L` fill the lower triangle of `LU`, while the elements of 
      * `U` fill `LU`'s upper triangle. The diagonal unity elements, which belong 
      * to `L`, are not stored.
      *
      * Note that this is NOT equal to the matrix `L * U` (which will be definition 
      * simply be a row-permutation of the original matrix `M`). 
      */
    SimpleMatrix<T,M,N> LU;
    
    /// The matrix `P` such that `P * M = L * U`, where `M` is the matrix decomposed within.
    PermutationMatrix<M> P;
    
protected:
    
    // LU stores both the upper and lower triangular parts of the decomposition.
    // The diagonal one elements (which belong to L) are not stored.
    bool singular;
    bool swap_parity;

public:
    
    /**
     * @brief Construct an empty PLU decomposition with dimensions `rows` x `cols`. 
     */
    PLUDecomposition(index_t rows, index_t cols):
            LU(rows, cols),
            P(rows),
            singular(false),
            swap_parity(false) {}
    /**
     * @brief Construct a PLU decompostion of `m`. `Mx` must be a matrix type.
     */
#ifdef PARSING_DOXYGEN
    template <typename Mx> explicit PLUDecompostion(const Mx &m) {}
#endif
    
    // what a garbage fire of a function signature :(
    template <typename Mx>
    explicit PLUDecomposition(
        const Mx& mx, 
        typename boost::enable_if_c
        <
            detail::MatrixDimensionMatch
            <
                SimpleMatrix<T,M,N>,
                Mx
            >::isStaticMatch,
            int
        >::type dummy=0):
            LU(mx.rows(), mx.cols()),
            P(mx.rows()),
            singular(false),
            swap_parity(false)
    {
        decompose(mx);
    }
    
    
    /**
     * @brief Decompose `mx` into this PLUDecomposition.
     *
     * Fill this PLU's buffers with the decomposition of `mx`.
     * 
     * `mx` must have the same number of rows and colums as `this.LU`.
     *
     * @param mx A matrix with dimensions matching `LU`. 
     * @return The number of degenerate rows in `mx`.
     */
    template <typename Mx>
    bool decompose(const Mx& mx) {
        // use P's existing buffer:
        index_t* P_i = P.getSrcData();
        std::copy(mx.begin(), mx.begin() + mx.rows() * mx.cols(), LU.begin());
        index_t degenerate_rows = decomp_plu(
            LU.begin(),
            LU.rows(),
            LU.cols(),
            P_i,
            &swap_parity);
        if (degenerate_rows > 0) singular = true;
        P.setRowSources(P_i);
        return degenerate_rows;
    }
    
    /**
     * @return The number of elements in the diagonal of `LU`. 
     */
    inline index_t diagonal() const {
        return std::min(LU.rows(), LU.cols());
    }
    
    /**
     * Get the lower triangular matrix.
     * 
     * @tparam S Element type of destination matrix.
     * @tparam J Row dimension of destination matrix. Must be 0 or equal to `LU.rows()`.
     * @tparam K Column dimension of destination matrix. Must be 0 or equal to `diagonal()`.
     * 
     * @param [out] into A matrix of dimension `LU.rows() x diagonal()`
     */
#ifdef PARSING_DOXYGEN
    template <typename S, index_t J, index_t K>
    void get_L(SimpleMatrix<S,J,K>* into) const {}
#endif
    template <typename S, index_t J, index_t K, MatrixLayout Lyt, StoragePolicy SP>
    inline typename boost::enable_if_c<
            detail::MatrixDimensionMatch<L_t, SimpleMatrix<S,J,K> >::isStaticMatch,
        void>::type 
    get_L(SimpleMatrix<S,J,K,Lyt,SP>* into) const {
#ifdef GEOMC_MTX_CHECK_DIMS
        const index_t diag = diagonal();
        if ((J * L_t::ROWDIM == 0 or J != L_t::ROWDIM or
             K * L_t::COLDIM == 0 or K != L_t::COLDIM) and 
             (into->rows() != LU.rows() or into->cols() != diag)) {
            throw DimensionMismatchException(into->rows(), into->cols(), LU.rows(), diag);
        }
        if ((J * K == 0 or J != K) and into->rows() != into->cols()) {
            throw NonsquareMatrixException(into->rows(), into->cols());
        }
#endif
        _copyL(into);
    }
    
    /**
     * Get the upper triangular matrix.
     * 
     * @tparam S Element type of destination matrix.
     * @tparam J Row dimension of destination matrix. Must be 0 or equal to `diagonal()`.
     * @tparam K Column dimension of destination matrix. Must be 0 or equal to `LU.cols()`.
     * 
     * @param [out] into A matrix of dimension `diagonal() x LU.cols()`
     */
#ifdef PARSING_DOXYGEN
    template <typename S, index_t J, index_t K>
    void get_U(SimpleMatrix<S,J,K> *into) const {}
#endif
    template <typename S, index_t J, index_t K, MatrixLayout Lyt, StoragePolicy SP>
    inline typename boost::enable_if_c<
            detail::MatrixDimensionMatch<U_t, SimpleMatrix<S,J,K> >::isStaticMatch,
        void>::type
    get_U(SimpleMatrix<S,J,K,Lyt,SP> *into) const {
#ifdef GEOMC_MTX_CHECK_DIMS
        const index_t diag = diagonal();
        if ((J * U_t::ROWDIM == 0 or J != U_t::ROWDIM or 
             K * U_t::COLDIM == 0 or K != U_t::COLDIM) and
             (into->rows() != diag or into->cols() != LU.cols())) {
            throw DimensionMismatchException(into->rows(), into->cols(), diag, LU.cols());
        }
#endif
        _copyU(into);
    }
    
    //EDIT: shouldn't it be U.cols()?
    // only vectors of dimension <U.rows()> are permitted.
    // we'll assume (perhaps unfairly?) that the client
    // is passing the right thing, since there is no
    // better way to enforce this.
    /**
     * Solve the linear matrix equation `Mx = b` for `x`,
     * where `M` is the square matrix decomposed herein. `b` must
     * be of length `LU.rows()`. 
     * @param [out] dest Array of length `LU.rows()`.
     * @param [in]  b Array of length `LU.rows()`.
     */
    template <typename S>
    inline void linear_solve(S *dest, const S *b) const {
        
#ifdef GEOMC_MTX_CHECK_DIMS
        _check_is_square();
#endif
        
#ifdef GEOMC_MTX_CHECK_ALIASING
        if (dest == b) {
            // because of the permutation, 
            // <b> will be destructively updated as it is read.
            index_t n = LU.rows();
            UniqueStorage<S,M> buf(n);
            std::copy(b, b+n, buf.get());
            geom::linear_solve_plu(
                LU.begin(),
                P.getRowSources(),
                LU.rows(),
                dest,
                buf.get());
            return;
        }
#endif
        
        geom::linear_solve_plu(LU.begin(), P.getRowSources(), LU.rows(), dest, b);
     }
    
    
    /**
     * Solve the linear matrix equation `Mx = b` for `x`, where `M` is the
     * square matrix decomposed herein.
     * 
     * @tparam S Element type of destination vector.
     * @tparam K Dimension of destination vector. Must be equal to `LU.rows()`.
     * 
     * @param b A vector of length `LU.rows()`.
     * @return The solution vector `x` to `Mx = b`.
     */
    template <typename S, index_t K>
    inline Vec<S,K> linear_solve(const Vec<S,K> &b) const {
        Vec<S,K> dest;

        // vectors get a dimension check, because we can:
#ifdef GEOMC_MTX_CHECK_DIMS
        if ((M == DYNAMIC_DIM or M != K) && LU.rows() != K) {
            throw DimensionMismatchException(LU.rows(), 1, K, 1);
        }
        _check_is_square();
#endif
        
        geom::linear_solve_plu(
            LU.begin(), 
            P.getRowSources(), 
            LU.rows(), 
            dest.begin(), 
            b.begin());
        
        return dest;
    }
    
    /**
     * Compute the matrix inverse of the decomposed matrix `M`, and copy it
     * to `into`. `M` must be square.
     *
     * @tparam J Row dimension of destination matrix. Must be 0 (dynamic) or `LU.rows()`.
     * @tparam K Column dimension of destination matrix. Must be 0 or `LU.cols()`.
     *  
     * @param [out] into Destination matrix; a square matrix with dimensions 
     * equal to `LU`.
     */
    template <index_t J, index_t K, MatrixLayout Lyt, StoragePolicy SP>
    void inverse(SimpleMatrix<T,J,K,Lyt,SP>* into) const {
        // todo: this can all probably be cleaner/cheaper
#ifdef GEOMC_MTX_CHECK_DIMS
        _check_is_square();
        // destination is valid?
        if ((J * K == 0 or J != M or K != N) and 
            (into->rows() != LU.rows() or into->cols() != LU.cols())) {
            throw DimensionMismatchException(
                into->rows(),
                into->cols(),
                LU.rows(),
                LU.cols());
        }
#endif
        T* dest = into->data_begin();
        
        // Set LUx = PI and solve for x, choosing columns of PI one at a time.
        // (If we transpose LU, also transpose PI)
        const index_t* p = P.getColSources();
        const index_t  n = LU.rows();
        std::fill(dest, dest + (n * n), 0);
        for (index_t i = 0; i < n; i++, dest += n) {
            dest[p[i]] = 1;
            geom::linear_solve_lu<T>(LU.begin(), LU.rows(), dest, dest);
        }
        if (Lyt == ROW_MAJOR) {
            // transpose the destination matrix
            auto& m = *into;
            for (index_t i = 0; i < n; i++) {
                for (index_t j = 0; j < i; j++) {
                    std::swap(m(i,j), m(j,i));
                }
            }
        }
    }
    
    
    /**
     * Compute the determinant of the decomposed matrix.
     */
    T det() const {
        const index_t n = LU.rows();
        T k = get_parity();
        for (index_t i = 0; i < n; i++) {
            k *= LU(i,i);
        }
        return k;
    }
    
    
    /**
     * @return `true` if the decomposed matrix is singular, `false` otherwise.
     */
    inline bool is_singular() const {
        return singular;
    }
    
    /**
     * @return -1 if `P` introduces a coordinate system handedness-swap; 1 otherwise.
     */
    inline int get_parity() const {
        return swap_parity ? -1 : 1;
    }
    
protected:
    
    template <typename Mx>
    void _copyL(Mx* into) const {
        for (index_t r = 0; r < LU.rows(); r++) {
            for (index_t c = 0; c < LU.cols(); c++) {
                T v;
                if (c < r) {
                    v = LU(r,c);
                } else if (c == r) {
                    v = 1;
                } else {
                    v = 0;
                }
                into->set(r,c,v);
            }
        }
    }
    
    template <typename Mx>
    void _copyU(Mx* into) const {
        for (index_t r = 0; r < LU.rows(); r++) {
            for (index_t c = 0; c < LU.cols(); c++) {
                T v;
                if (c >= r) {
                    v = LU(r,c);
                } else {
                    v = 0;
                }
                into->set(r,c,v);
            }
        }
    }
    
    inline void _check_is_square() const {
        if ((M * N == 0 or M != N) and LU.rows() != LU.cols()) {
            throw NonsquareMatrixException(LU.rows(), LU.cols());
        }
    }
    
};  // plu class

}  // namespace geom

#endif /* LUDECOMP_H_ */
