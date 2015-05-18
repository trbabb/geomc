/*
 * DiagMatrix.h
 *
 *  Created on: Apr 28, 2013
 *      Author: tbabb
 */

#ifndef DIAGMATRIX_H_
#define DIAGMATRIX_H_

#include <geomc/linalg/mtxdetail/MatrixBase.h>

namespace geom {

namespace detail {

    // const T because m[r][c].mutate() does not make sense. m[r][c] = X is allowed, however.
    template <typename T, index_t M, index_t N>
    struct _ImplMtxReftype<DiagMatrix<T,M,N>,T> {
        typedef detail::MtxAssignmentProxy<DiagMatrix<T,M,N>, const T> reference;
    };
    
}; // end namespace detail

/** @ingroup matrix 
 *  @brief A matrix with nonzero elements only along the main diagonal.
 * 
 * @tparam T Element type.
 * @tparam M Row dimension.
 * @tparam N Column dimension.
 * 
 * Storage for `DiagMatrix`es is O(n), and operations with diagonal matrices 
 * usually benefit from the O(n) algorithmic improvement associated with their 
 * lower dimension.
 */
template <typename T, index_t M, index_t N>
class DiagMatrix : public detail::MatrixBase<T,M,N, DiagMatrix<T,M,N> > {
    
    index_t diag;
#ifndef PARSING_DOXYGEN
    // this shit confuses doxygen because it am stoopid
    UniqueStorage<T, ((M<N)?M:N) > data;
#endif
    typename Dimension<M>::storage_t n_rows;
    typename Dimension<N>::storage_t n_cols;
    
public:
    
    typedef typename detail::_ImplMtxReftype<DiagMatrix<T,M,N>,T>::reference reference;
    
    using detail::MatrixBase<T,M,N, DiagMatrix<T,M,N> >::get;
    
    /**
     * Construct a new diagonal matrix.
     * @param nrows Number of rows (ignored / not required if statically sized).
     * @param ncols Number of columns (ignored / not required if statically sized).
     */
#ifdef PARSING_DOXYGEN
    explicit DiagMatrix(index_t nrows, index_t ncols) {}
#else
    explicit DiagMatrix(index_t nrows=detail::DefinedIf<M != DYNAMIC_DIM, M>::value, 
               index_t ncols=detail::DefinedIf<N != DYNAMIC_DIM, N>::value) :
                     diag(nrows < ncols ? nrows : ncols),
                     data(diag) {
        Dimension<M>::set(n_rows, nrows);
        Dimension<N>::set(n_cols, ncols);
        setIdentity(); 
    }
#endif
    
    /**
     * Construct a new `n` x `n` diagonal matrix, with diagonal elements copied 
     * from `src`.
     * @param src Array containing diagonal elements.
     * @param n Number of elements in `src`. Ignored / not required if statically sized.
     */
#ifdef PARSING_DOXYGEN
    DiagMatrix(const T src[], index_t n) {}
#else
    DiagMatrix(const T src[], 
               index_t n=detail::DefinedIf<M * N != DYNAMIC_DIM, (M < N ? M : N)>::value) :
                   diag(n) {
        Dimension<M>::set(n_rows, n);
        Dimension<N>::set(n_cols, n);
        std::copy(src, src + n, data.get());
    }
#endif
    
    inline index_t rows() const {
        return Dimension<M>::get(n_rows);
    }
    
    inline index_t cols() const {
        return Dimension<N>::get(n_cols);
    }
    
    inline T get(index_t r, index_t c) const {
#ifdef GEOMC_MTX_CHECK_BOUNDS
            if (r < 0 or r >= rows() or c < 0 or c >= cols()) {
                throw std::out_of_range("matrix coordinates");
            }
#endif
        return (r == c) ? data.get()[r] : (T(0));
    }

    void setIdentity() {
        std::fill(data.get(), data.get() + diag, 1);
    }
    
    inline void setZero() {
        std::fill(data.get(), data.get() + diag, 0);
    }
    
    /**
     * @return An iterator over the elements of the main diagonal (those with
     * matching row and column coordinates), pointing at element `(0, 0)`. 
     */
    T* diagonal_begin() {
        return data.get();
    }
    
    /**
     * @return A read-only iterator over the elements of the main diagonal (those with
     * matching row and column coordinates), pointing at element `(0, 0)`. 
     */
    const T* diagonal_begin() const {
        return data.get();
    }
    
    /**
     * @return An iterator over the elements of the main diagonal (those with
     * matching row and column coordinates), pointing at the last diagonal element.
     */
    T* diagonal_end() {
        return data.get() + diag;
    }
    
    /**
     * @return A read-only iterator over the elements of the main diagonal (those with
     * matching row and column coordinates), pointing at the last diagonal element.
     */
    const T* diagonal_end() const {
        return data.get() + diag;
    }

    /**
     * Return the `i`th diagonal element.
     */
    const T& diagonal(index_t i) const {
        return data.get()[i];
    }

    /**
     * Return a reference to the `i`th diagonal element.
     */
    T& diagonal(index_t i) {
        return data.get()[i];
    }
    
    inline void getStorageIDs(storage_id_t *buf) const {
        *buf = data.get();
    }
    
    inline index_t getStorageIDCount() const {
        return 1;
    }
};

}; // end namespace geom

#endif /* DIAGMATRIX_H_ */
