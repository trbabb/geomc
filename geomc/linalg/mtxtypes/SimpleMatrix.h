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
#include <geomc/linalg/mtxdetail/MatrixCopy.h>

// TODO: templatize to intermediate class, and make a child class
//       which is a wrapper for user-owned memory.

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

/** @ingroup matrix 
 *  @brief A basic matrix with `M x N` elements. 
 * 
 * @tparam T Element type.
 * @tparam M Row dimension.
 * @tparam N Column dimension.
 * 
 * If `M` or `N` are 0, the matrix has runtime-chosen size, and all copy-constructed
 * duplicates of this matrix should be treated as references to a common array.
 * 
 *  For more on matrices, see the @link matrix matrix module documentation@endlink.
 */
template <typename T, index_t M, index_t N>
class SimpleMatrix : public detail::WriteableMatrixBase<T,M,N, SimpleMatrix<T,M,N> > {
    Storage<T,M*N> data;
    typename Dimension<M>::storage_t n_rows; // this contortion prevents allocating storage 
    typename Dimension<N>::storage_t n_cols; // for a dimension size when it is not dynamic.
    
public:
    
    // MatrixBase types
    /// @brief Read-only row-major iterator over the elements of the matrix.
    typedef const T* const_iterator;
    /// @brief Read-only iterator over the elements of a row.
    typedef const T* const_row_iterator;
    
    // WriteableMatrixBase types
    /// @brief Reference to element type.
    typedef T& reference;
    /// @brief Writeable row-major iterator over matrix elements.
    typedef T* iterator;
    /// @brief Writeable iterator over elements of a row.
    typedef T* row_iterator;
    
    /**
     * @brief Construct a matrix of size `(nrows x ncols)`. 
     * 
     * For matrices with dynamic dimension, a size argument is **required** for 
     * that dimension. Size arguments will be ignored for dimensions that are 
     * statically-sized.
     * 
     * Examples:
     *     
     *     SimpleMatrix<double, 3, 3> m1;
     *     SimpleMatrix<double, 0, 0> m2(3, 3);
     *     SimpleMatrix<double, 0, 0> m3; // XXX: Compiler error!
     * 
     * @param nrows Number of rows in the matrix.
     * @param ncols Number of columns in the matrix.
     */
#ifdef PARSING_DOXYGEN
    explicit SimpleMatrix(index_t nrows, index_t ncols) {}
#else
    // I would not expect this trickery to work, but it does!
    // this mandates row, column arguments for dynamic matrices.
    // arguments are basically ignored if matrix is static 
    explicit SimpleMatrix(index_t nrows=detail::DefinedIf<M != DYNAMIC_DIM, M>::value, 
                 index_t ncols=detail::DefinedIf<N != DYNAMIC_DIM, N>::value) : 
                     data(nrows * ncols) {
        Dimension<M>::set(n_rows, nrows);
        Dimension<N>::set(n_cols, ncols);
        setIdentity(); 
    }
#endif

    /**
     * @brief Construct a matrix of size `(nrows x ncols)`, initialized with `src_data`. 
     * 
     * For matrices with dynamic dimension, a size argument is **required** for 
     * that dimension. Size arguments will be ignored for dimensions that are 
     * statically-sized.
     * 
     * @param src_data Array of `nrows * ncols` elements, in row-major order, to
     * be copied.
     * @param nrows Number of rows in the matrix.
     * @param ncols Number of columns in the matrix.
     */
#ifdef PARSING_DOXYGEN
    explicit SimpleMatrix(const T* src_data, index_t nrows, index_t ncols) {}
#else
    explicit SimpleMatrix(
                 const T* src_data,
                 index_t nrows=detail::DefinedIf<M != DYNAMIC_DIM, M>::value, 
                 index_t ncols=detail::DefinedIf<N != DYNAMIC_DIM, N>::value) : 
                     data(nrows * ncols) {
        Dimension<M>::set(n_rows, nrows);
        Dimension<N>::set(n_cols, ncols);
        std::copy(src_data, src_data + nrows * ncols, data.get());
    }
#endif
    
    /**
     * Construct and initialize this matrix with the contents of another.
     * @tparam Mx A matrix type with agreeing dimension. 
     * @param mtx Matrix containing source elements.
     */
#ifdef PARSING_DOXYGEN
    template <typename Mx>
    SimpleMatrix(const Mx &mtx) {}
#else
    template <typename Mx>
    SimpleMatrix(const Mx &mtx,
                 typename boost::enable_if_c<
                    (detail::LinalgDimensionMatch<SimpleMatrix<T,M,N>, Mx>::val),
                    void*>::type dummy=0) {
        Dimension<M>::set(n_rows, mtx.rows());
        Dimension<N>::set(n_cols, mtx.cols());
        detail::_mtxcopy(this, mtx);
    }
#endif
    
    /**************************
     * Operators              *
     **************************/
    
    /**
     * Scalar in-place multiplication.
     * @tparam U  A type satisfying `boost::is_scalar<>`.
     * @param  k  Scalar value.
     * @return    A reference to `this`. 
     */
#ifdef PARSING_DOXYGEN
    template <typename U> SimpleMatrix<T,M,N>& operator*=(U k) {}
#else
    template <typename U>
    typename boost::enable_if<boost::is_scalar<U>, SimpleMatrix<T,M,N>&>::type
    operator*=(U k) {
        T *p = data.get();
        for (index_t i = 0; i < rows() * cols(); i++) {
            p[i] *= k;
        }
        return *this;
    }
#endif
    
    /**************************
     * Functions              *
     **************************/
    
    /**
     * @return Number of rows in the matrix.
     */
    inline index_t rows() const {
        return Dimension<M>::get(n_rows);
    }
    
    /**
     * @return number of columns in the matrix.
     */
    inline index_t cols() const {
        return Dimension<N>::get(n_cols);
    }
    
    /**
     * @param i Index of row (zero-indexed)
     * @return A read-only iterator over the elements of row `i`.
     */
    inline const_row_iterator row(index_t i) const {
        return data.get() + (cols() * i);
    }
    
    /**
     * @param i Index of row (zero-indexed)
     * @return A writeable iterator over the elements of row `i`.
     */
    inline row_iterator row(index_t i) {
        return data.get() + (cols() * i);
    }
    
    /**
     * @return A read-only, random-access, row-major iterator over the elements of this matrix,
     * pointing to the element at (0,0).
     */
    inline const_iterator begin() const {
        return data.get();
    }
    
    /**
     * @return A read-only, random-access, row-major iterator over the elements of this matrix,
     * pointing to the just beyond the last element in the matrix.
     */
    inline const_iterator end() const {
        return data.get() + (rows() * cols());
    }
    
    /**
     * @return A writeable, random-access, row-major iterator over the elements of this matrix,
     * pointing to the element at (0,0).
     */
    inline iterator begin() {
        return data.get();
    }
    
    /**
     * @return A writeable, random-access, row-major iterator over the elements of this matrix,
     * pointing to the element just beyond the last element in the matrix.
     */
    inline iterator end() {
        return data.get() + (rows() * cols());
    }
    
    /**
     * Get the element at `(row, col)`.
     * @param r Zero-indexed row coordinate
     * @param c Zero-indexed column coordinate
     * @return The element at `(row, col)`
     */
    inline T get(index_t r, index_t c) const {
        return data.get()[r * cols() + c];
    }
    
    /**
     * Get the element at `(row, col)`.
     * @param r Zero-indexed row coordinate
     * @param c Zero-indexed column coordinate
     * @return A reference to the element at `(row, col)`
     */
    inline reference get(index_t r, index_t c) {
        return data.get()[r * cols() + c];
    }
    
    /**
     * Set the element at `(row, col)` to `val`.
     * @param r Zero-indexed row coordinate
     * @param c Zero-indexed column coordinate
     * @param val New value of element at `(row, col)`
     * @return A reference to the element at `(row, col)`, for convenience.
     */
    inline reference set(index_t r, index_t c, T val) {
        return (data.get()[r * cols() + c] = val);
    }
    
    /**
     * In-place transpose this matrix.
     */
    void transpose() {
        for (int r = 0; r < rows(); r++) {
            for (int c = 0; c < r; c++) {
                std::swap(get(r,c), get(c,r));
            }
        }
    }

    /**
     * Clear this matrix and set its elements to the identity matrix.
     */
    void setIdentity() {
        setZero();
        T *last = end();
        for (T *p = begin(); p < last; p += cols() + 1) {
            *p = 1;
        }
    }
    
    /**
     * Set all the elements of this matrix to zero.
     */
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
