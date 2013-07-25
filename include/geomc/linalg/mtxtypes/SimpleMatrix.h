/*
 * SimpleMatrix.h
 *
 *  Created on: Apr 18, 2013
 *      Author: tbabb
 */

#ifndef SIMPLEMATRIX_H_
#define SIMPLEMATRIX_H_

#include <boost/utility/enable_if.hpp>
#include <geomc/linalg/mtxdetail/MatrixBase.h>

// TODO: templatize across row/col-major layout.
//       ...or just make a whole new class.

namespace geom {

namespace detail {

    // specializations of the iterator type helper classes
    // tell the base class which iterators to use.

    template <typename T, index_t M, index_t N>
    struct _ImplMtxReftype<SimpleMatrix<T,M,N>,T> {
        typedef T& reference;
    };
    
    template <typename T, index_t M, index_t N>
    struct _ImplMtxRowtype<SimpleMatrix<T,M,N>,T> {
        typedef T* row_iterator;
    };
    
    template <typename T, index_t M, index_t N>
    struct _ImplMtxConstRowtype<SimpleMatrix<T,M,N>,T> {
        typedef const T* const_row_iterator;
    };
    
}; // end namespace detail


template <typename T, index_t M, index_t N>
class SimpleMatrix : public detail::WriteableMatrixBase<T,M,N, SimpleMatrix<T,M,N> > {
    detail::Storage<T,M*N> data;
    typename detail::Dimension<M>::storage_t n_rows; // this contortion prevents allocating storage 
    typename detail::Dimension<N>::storage_t n_cols; // for a dimension size when it is not dynamic.
    
public:
    
    // MatrixBase types
    typedef const T* const_iterator;
    typedef const T* const_row_iterator;
    
    // WriteableMatrixBase types
    typedef T& reference;
    typedef T* iterator;
    typedef T* row_iterator;
    
    // I would not expect this trickery to work, but it does!
    // this mandates row, column arguments for dynamic matrices.
    // arguments are basically ignored if matrix is static 
    SimpleMatrix(index_t nrows=detail::DefinedIf<M != DYNAMIC_DIM, M>::value, 
                 index_t ncols=detail::DefinedIf<N != DYNAMIC_DIM, N>::value) : 
                     data(nrows * ncols) {
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
    
    inline const_row_iterator row(index_t i) const {
        return data.get() + (cols() * i);
    }
    
    inline row_iterator row(index_t i) {
        return data.get() + (cols() * i);
    }
    
    inline const_iterator begin() const {
        return data.get();
    }
    
    inline const_iterator end() const {
        return data.get() + (rows() * cols());
    }
    
    inline iterator begin() {
        return data.get();
    }
    
    inline iterator end() {
        return data.get() + (rows() * cols());
    }
    
    inline T get(index_t r, index_t c) const {
        return data.get()[r * cols() + c];
    }
    
    inline reference get(index_t r, index_t c) {
        return data.get()[r * cols() + c];
    }
    
    inline reference set(index_t r, index_t c, T val) {
        return (data.get()[r * cols() + c] = val);
    }
    
    void transpose() {
        for (int r = 0; r < rows(); r++) {
            for (int c = 0; c < r; c++) {
                std::swap(get(r,c), get(c,r));
            }
        }
    }

    void setIdentity() {
        setZero();
        T *last = end();
        for (T *p = begin(); p < last; p += cols() + 1) {
            *p = 1;
        }
    }
    
    inline void setZero() {
        std::fill(begin(), end(), 0);
    }
    
    inline void getStorageIDs(storage_id_t *buf) const {
        *buf = data.get();
    }
    
    inline index_t getStorageIDCount() const {
        return 1;
    }
};


}; // end namespace geom

#endif /* SIMPLEMATRIX_H_ */
