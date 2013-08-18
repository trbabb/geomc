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
 * @brief A matrix with nonzero elements only along the main diagonal.
 * 
 * Storage for `DiagMatrix`es is O(n), and operations with diagonal matrices 
 * usually benefit from the O(n) algorithmic improvement associated with their 
 * lower dimension.
 */
template <typename T, index_t M, index_t N>
class DiagMatrix : public detail::WriteableMatrixBase<T,M,N, DiagMatrix<T,M,N> > {
    index_t diag;
    detail::Storage<T, ((M<N)?M:N) > data;
    typename detail::Dimension<M>::storage_t n_rows;
    typename detail::Dimension<N>::storage_t n_cols;
    
public:
    
    typedef typename detail::_ImplMtxReftype<DiagMatrix<T,M,N>,T>::reference reference;
    
    using detail::WriteableMatrixBase<T,M,N, DiagMatrix<T,M,N> >::get;
    
    DiagMatrix(index_t nrows=detail::DefinedIf<M != DYNAMIC_DIM, M>::value, 
               index_t ncols=detail::DefinedIf<N != DYNAMIC_DIM, N>::value) :
                     diag(nrows < ncols ? nrows : ncols),
                     data(diag) {
        detail::Dimension<M>::set(n_rows, nrows);
        detail::Dimension<N>::set(n_cols, ncols);
        setIdentity(); 
    }
    
    inline index_t rows() const {
        return detail::Dimension<M>::value(n_rows);
    }
    
    inline index_t cols() const {
        return detail::Dimension<N>::value(n_cols);
    }
    
    inline T get(index_t r, index_t c) const {
        #ifdef GEOMC_MTX_CHECK_BOUNDS
            if (r < 0 or r >= rows() or c < 0 or c >= cols()) {
                throw std::out_of_range("matrix coordinates");
            }
        #endif
        return (r == c) ? data.get()[r] : (T(0));
    }
    
    inline reference set(index_t r, index_t c, T val) {
        if (r == c && c < diag) {
            data.get()[r] = val;
        } else if (val != 0) {
            throw std::out_of_range("off-diagonal write");
        }
        return reference(this, r, c);
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
    
    inline void getStorageIDs(storage_id_t *buf) const {
        *buf = data.get();
    }
    
    inline index_t getStorageIDCount() const {
        return 1;
    }
};

}; // end namespace geom

#endif /* DIAGMATRIX_H_ */
