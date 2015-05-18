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

/************************************************
 * Matrix LU decomposition                      *
 *                                              *
 * Implementation using Doolittle's algorithm   *
 ************************************************/

#define _MxElem(r,c) m[cols*r + c]
    
/**
 * LU decomposition, pivoting columns.
 * 
 * @param m Row-major array of elements to decompose.
 * @param rows Number of rows in `m`.
 * @param cols Number of columns in `m`.
 * @param reorder Array with space for `cols` elements to be filled with the column source indexes.
 * @param swap_parity Whether an odd number of column-swaps was performed.
 * @return The number of degenerate columns discovered.
 * @related PLUDecomposition
 */
template <typename T>
index_t decompLUP(T* m, index_t rows, index_t cols, index_t *reorder, bool *swap_parity) {
    const index_t n = std::min(rows, cols);
    index_t degenerate_ct = 0;
    *swap_parity = false;
    // fill permutation array
    for (index_t i = 0; i < cols; i++){
        reorder[i] = i;
    }
    
    for (index_t i = 0; i < n - 1; i++) {
        // find pivot
        T biggest = std::abs(_MxElem(i,i));
        index_t pvt = i;
        for (index_t p = i + 1; p < cols; p++) {
            T pvt_val = std::abs(_MxElem(i,p));
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
                std::swap(_MxElem(r, i),
                          _MxElem(r, pvt));
            }
            // make a note of the permutation in <P>
            std::swap(reorder[i], reorder[pvt]);
            *swap_parity = !*swap_parity;
        }
        
        // eliminate lower elements
        T a = _MxElem(i,i);
        for (index_t r = i + 1; r < rows; r++) {
            T b = _MxElem(r,i) / a;
            for (index_t c = i + 1; c < cols; c++) {
                // R_r = R_r - b * R_i
                T src_elem = _MxElem(i,c);
                T dst_elem = _MxElem(r,c);
                _MxElem(r,c) = dst_elem - b * src_elem;
            }
            // set the lower matrix
            _MxElem(r,i) = b;
        }
    }
    
    return degenerate_ct;
}

// we operate on a bare data array (in row-contiguous order).
// this saves on multiple instantiations of this method for
// different static sizes of SimpleMatrix. It's 7% faster too!
/**
 * LU decomposition, pivoting rows.
 * 
 * @param m Row-major array of elements to decompose.
 * @param rows Number of rows in `m`.
 * @param cols Number of columns in `m`.
 * @param reorder Array with space for `rows` elements to be filled with the row source indexes.
 * @param swap_parity Whether an odd number of row-swaps was performed.
 * @return The number of degenerate rows discovered.
 * @related PLUDecomposition
 */
template <typename T>
bool decompPLU(T* m, index_t rows, index_t cols, index_t *reorder, bool *swap_parity) {
    const index_t n = std::min(rows, cols);
    index_t degenerate_ct = 0;
    *swap_parity = false;
    // fill permutation array
    for (index_t i = 0; i < rows; i++){
        reorder[i] = i;
    }
    
    for (index_t i = 0; i < n - 1; i++) {
        // find pivot
        T biggest = std::abs(_MxElem(i,i));
        index_t pvt = i;
        for (index_t p = i + 1; p < rows; p++) {
            T pvt_val = std::abs(_MxElem(p,i));
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
            // swap row <i> with row <pvt>
            for (index_t c = 0; c < cols; c++) {
                std::swap(_MxElem(i,   c),
                          _MxElem(pvt, c));
            }
            // make a note of the permutation in <P>
            std::swap(reorder[i], reorder[pvt]);
            *swap_parity = !*swap_parity;
        }
        
        // eliminate lower elements
        T a = _MxElem(i,i);
        for (index_t r = i + 1; r < rows; r++) {
            T b = _MxElem(r,i) / a;
            for (index_t c = i + 1; c < cols; c++) {
                // R_r = R_r - b * R_i
                T src_elem = _MxElem(i,c);
                T dst_elem = _MxElem(r,c);
                _MxElem(r,c) = dst_elem - b * src_elem;
            }
            // set the lower matrix
            _MxElem(r,i) = b;
        }
    }
    
    return degenerate_ct;
}

#undef _MxElem

//////////// PLU class ////////////

/** 
 * @ingroup matrix
 * @brief Computes the PLU decompostion for a matrix `A`, such that `PA = LU`.
 * 
 * `L` and `U` are lower and upper triangular matrices, respectively.
 * `P` has dimension `(LU.rows() x LU.rows())`, and `LU` has the dimension of `A`.
 */
template <typename T, index_t M, index_t N>
class PLUDecomposition {
public:
    /// Dimension of the diagonal of `LU`. The minimum of `M` and `N`, or 0 if either dimension is dynamic.
    static const index_t DIAG = M<N?M:N;
    /// Matrix type for holding `L`.
    typedef SimpleMatrix<T,M,DIAG> L_t;
    /// Matrix type for holding `U`.
    typedef SimpleMatrix<T,DIAG,N> U_t;
    
protected:
    
    // LU stores both the upper and lower triangular parts of the decomposition.
    // The diagonal one elements (which belong to L) are not stored.

    SimpleMatrix<T,M,N> LU;
    PermutationMatrix<M> P;
    bool singular;
    bool swap_parity;
    
    PLUDecomposition(index_t n_r, index_t n_c):
            LU(n_r, n_c),
            P(n_r),
            singular(false),
            swap_parity(false) {}
    
public:
    /**
     * Construct a PLU decompostion of `m`. `Mx` must be a matrix type.
     */
#ifdef PARSING_DOXYGEN
    template <typename Mx> explicit PLUDecompostion(const Mx &m){}
#endif
    
    template <typename Mx>
    explicit PLUDecomposition(const Mx& m, 
                              typename boost::enable_if_c<detail::MatrixDimensionMatch<SimpleMatrix<T,M,N>, Mx>::isStaticMatch, int>::type dummy=0):
            LU(m.rows(), m.cols()),
            P(m.rows()),
            singular(false),
            swap_parity(false) {
        
        mtxcopy(&LU, m);
        // TODO: this alloc isn't necessary if we become a friend of PermutationMatrix
        UniqueStorage<index_t, Mx::ROWDIM> reorder(m.rows());
        bool ok = decompPLU(LU.begin(), m.rows(), m.cols(), reorder.get(), &swap_parity) == 0;

        if (not ok) {
            singular = true;
        }
        P.setRowSources(reorder.get());
    }
    
public:
    
    /**
     * @return The number of elements in the diagonal of `LU`. 
     */
    inline index_t diagonal() const {
        return std::min(LU.rows(), LU.cols());
    }
    
    /**
     * @return The row-permutation matrix.
     */
    const PermutationMatrix<M>& getP() const {
        return P;
    }
    
    /**
     * Get `L` and `U` as superimposed matrices. The elements of `L` fill
     * the lower triangle, and the elements of `U` fill the upper. The diagonal
     * unity elements, which belong to `L`, are not stored.
     */
    const SimpleMatrix<T,M,N>& getLU() const {
        return LU;
    }
    
    /**
     * @return A copy of the lower-triangular matrix.
     */
    const SimpleMatrix<T,M,DIAG> getL() const {
        typedef detail::_ImplMtxInstance< SimpleMatrix<T,M,DIAG> > instancer;
        const index_t diag = diagonal();
        SimpleMatrix<T,M,DIAG> out = instancer::instance(LU.rows(), diag);
        _copyL(&out);
        return out;
    }
    
    /**
     * @return A copy of the upper triangular matrix.
     */
    const SimpleMatrix<T,DIAG,N> getU() const {
        typedef detail::_ImplMtxInstance< SimpleMatrix<T,DIAG,N> > instancer;
        const index_t diag = diagonal();
        SimpleMatrix<T,DIAG,N> out = instancer::instance(diag, LU.cols());
        _copyU(&out);
        return out;
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
    void getL(SimpleMatrix<S,J,K> *into) const {}
#endif
    template <typename S, index_t J, index_t K, StoragePolicy SP>
    inline typename boost::enable_if_c<
            detail::MatrixDimensionMatch<L_t, SimpleMatrix<S,J,K> >::isStaticMatch,
        void>::type 
    getL(SimpleMatrix<S,J,K,SP> *into) const {
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
    void getU(SimpleMatrix<S,J,K> *into) const {}
#endif
    template <typename S, index_t J, index_t K, StoragePolicy SP>
    inline typename boost::enable_if_c<
            detail::MatrixDimensionMatch<U_t, SimpleMatrix<S,J,K> >::isStaticMatch,
        void>::type
    getU(SimpleMatrix<S,J,K,SP> *into) const {
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
    inline void linearSolve(S *dest, const S *b) const {
        
#ifdef GEOMC_MTX_CHECK_DIMS
        _checkIsSquare();
#endif
        
#ifdef GEOMC_MTX_CHECK_ALIASING
        if (dest == b) {
            // because of the permutation, <b> will be destructively
            // updated as it is read.
            index_t n = LU.rows();
            UniqueStorage<S,M> buf(n);
            std::copy(b, b+n, buf.get());
            _linearSolve(dest, buf.get());
            return;
        }
#endif
        
        _linearSolve(dest, b);
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
    inline Vec<S,K> linearSolve(const Vec<S,K> &b) const {
        Vec<S,K> dest;

        // vectors get a dimension check, because we can:
#ifdef GEOMC_MTX_CHECK_DIMS
        if ((M == DYNAMIC_DIM or M != K) && LU.rows() != K) {
            throw DimensionMismatchException(LU.rows(), 1, K, 1);
        }
        _checkIsSquare();
#endif
        
        _linearSolve(dest.begin(), b.begin());
        return dest;
    }
    
    /**
     * Compute the matrix inverse of the decomposed matrix `M`, and copy it
     * to `into`. `M` must be square.
     * 
     * @tparam S Element type of destination matrix.
     * @tparam J Row dimension of destination matrix. Must be 0 (dynamic) or `LU.rows()`.
     * @tparam K Column dimension of destination matrix. Must be 0 or `LU.cols()`.
     *  
     * @param [out] into Destination matrix; a square matrix with dimensions 
     * equal to `LU`.
     */
    template <typename S, index_t J, index_t K, StoragePolicy SP>
    void inverse(SimpleMatrix<S,J,K,SP> *into) const {
        
#ifdef GEOMC_MTX_CHECK_DIMS
        _checkIsSquare();
        // destination is valid?
        if ((J*K == 0 or J != M or K != N) and 
            (into->rows() != LU.rows() or into->cols() != LU.cols())) {
            throw DimensionMismatchException(into->rows(), into->cols(), LU.rows(), LU.cols());
        }
#endif
        _inverseTranspose(into->begin());
        into->transpose();
    }
    
    /**
     * Compute the determinant of the decomposed matrix.
     */
    T det() const {
        const index_t n = LU.rows();
        T k = getParity();
        for (index_t i = 0; i < n; i++) {
            k *= LU.get(i,i);
        }
        return k;
    }
    
    /**
     * @return `true` if the decomposed matrix is singular, `false` otherwise.
     */
    inline bool isSingular() const {
        return singular;
    }
    
    /**
     * @return -1 if `P` introduces a coordinate system handedness-swap; 1 otherwise.
     */
    inline int getParity() const {
        return swap_parity ? -1 : 1;
    }
    
protected:
    
    template <typename Mx>
    void _copyL(Mx *into) const {
        for (index_t r = 0; r < LU.rows(); r++) {
            for (index_t c = 0; c < LU.cols(); c++) {
                T v;
                if (c < r) {
                    v = LU.get(r,c);
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
    void _copyU(Mx *into) const {
        for (index_t r = 0; r < LU.rows(); r++) {
            for (index_t c = 0; c < LU.cols(); c++) {
                T v;
                if (c >= r) {
                    v = LU.get(r,c);
                } else {
                    v = 0;
                }
                into->set(r,c,v);
            }
        }
    }
    
    inline void _checkIsSquare() const {
        if ((M * N == 0 or M != N) and LU.rows() != LU.cols()) {
            throw NonsquareMatrixException(LU.rows(), LU.cols());
        }
    }
    
    
    template <typename S>
    void _linearSolve(S *dest, const S *b, bool permute=true) const {
        const index_t  n = LU.rows();
        const index_t *p = P.getRowSources(); // permutation map
        
        // <y> and <dest>'s elements are used such that
        // y[i] is never read after dest[i] is written.
        // thus we may collapse their storage and save space:
        S *y = dest;

        // LUx = Pb
        // (Ux) is a vector, so let's solve for it:
        // Ly  = Pb
        y[0] = permute ? b[p[0]] : b[0];
        for (index_t r = 1; r < n; r++) {
            y[r] = permute ? b[p[r]] : b[r];
            for (index_t c = 0; c < r; c++) {
                y[r] -= y[c] * LU.get(r,c);
            }
        }
        // now with y, we may obtain x from:
        // Ux = y
        for (index_t r = n - 1; r >= 0; r--) {
            // dest[r] = y[r]; // nop; dest and y are the same!
            for (index_t c = n - 1; c > r; c--) {
                dest[r] -= dest[c] * LU.get(r,c);
            }
            dest[r] /= LU.get(r,r);
        }
    }
    
    
    template <typename S>
    void _inverseTranspose(S *dest) const {
        // Set LUx = PI and solve for x, choosing columns of PI one at a time.
        // Here, we use the rows of our destination matrix as though they are
        // column vectors, and the caller will transpose.
        const index_t *p = P.getColSources();
        const index_t n = LU.rows();
        std::fill(dest, dest + (n*n), 0);
        for (index_t i = 0; i < n; i++, dest += n) {
            dest[p[i]] = 1;
            _linearSolve(dest, dest, false);
        }
    }
    
};  // plu class

}  // namespace geom

#endif /* LUDECOMP_H_ */
