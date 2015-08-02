/*
 * SimpleMatrix.h
 *
 * This involves a gambit very similar to the VecBase scenario:
 * We have a class which is essentially unchanging in its body, but
 * needs to have permuted constructors based on a template parameter.
 * Therefore we put the entire body in a base class, and then make a thin
 * derived class with partial specialization, and design the constructor
 * interface as we please for each specialization.
 *
 * This class structure is pretty darn ugly, so we hide it from the user
 * for the purposes of Doxygen by simply pretending the (hidden) common base class 
 * is the user-facing one.
 *
 *  Created on: Apr 18, 2013
 *      Author: tbabb
 */

#ifndef SIMPLEMATRIX_H_
#define SIMPLEMATRIX_H_

#include <boost/utility/enable_if.hpp>
#include <geomc/linalg/mtxdetail/MatrixBase.h>
#include <geomc/linalg/mtxdetail/MatrixCopy.h>


// TODO: templatize across row/col-major layout.
//       ...or just make a whole new class.

namespace geom {


///////////////////// Details /////////////////////


namespace detail {


    // allow myself to introduce myself:
    // I am a curiously recurring class, and I have to declare myself to declare myself.
    template <typename T, index_t M, index_t N, StoragePolicy P>
    class FlatMatrixBase;


    // specializations of the iterator type helper classes
    // tell the base class which iterators to use.

    template <typename T, index_t M, index_t N, StoragePolicy P>
    struct _ImplMtxReftype<FlatMatrixBase<T,M,N,P>,T> {
        typedef T& reference;
    };
    
    template <typename T, index_t M, index_t N, StoragePolicy P>
    struct _ImplMtxRowtype<FlatMatrixBase<T,M,N,P>,T> {
        typedef T* row_iterator;
    };
    
    template <typename T, index_t M, index_t N, StoragePolicy P>
    struct _ImplMtxConstRowtype<FlatMatrixBase<T,M,N,P>,T> {
        typedef const T* const_row_iterator;
    };
    
}; // end namespace detail



///////////////////// Flat matrix base class /////////////////////


#ifndef PARSING_DOXYGEN
namespace detail {
#endif


/** @ingroup matrix 
 *  @brief A basic matrix with `M x N` elements. 
 * 
 * @tparam T Element type.
 * @tparam M Row dimension.
 * @tparam N Column dimension.
 * @tparam P Policy for memory ownership of the backing array.
 *
 * Example:
 *
 *     SimpleMatrix<double,4,4> mx;
 * 
 * If `M` or `N` are 0, the matrix has runtime-chosen size.
 *
 * The storage policy behavior is as follows:
 * <ul>
 * <li>If the StoragePolicy is `STORAGE_UNIQUE` (default), then copy-constructed 
 * duplicates of dynamically-sized matrices will make a full copy of the underlying array.</li>
 * <li>If the StoragePolicy is `STORAGE_SHARED`, then all copy-constructed
 * duplicates of dynamically-sized matrixes should be treated as references to a common array.</li>
 * <li>If the StoragePolicy is `STORAGE_USER_OWNED`, then the matrix will wrap a user-provided
 * backing array, whose lifetime is managed manually.</li>
 * </ul>
 *
 * Note that in c++11, array duplicatations will use rvalue references to avoid a performing full 
 * array-copy where possible.
 *
 * Wrapping an array
 * -----------------
 *
 * A SimpleMatrix can be made into a wrapper around a user-owned array like so:
 *
 *     double myRowMajorArray[16] = { ... };
 *     SimpleMatrix<double,4,4,STORAGE_USER_OWNED> mtx(myRowMajorArray);
 *
 * In the above example, no array duplication will occur, and `myRowMajorArray` is 
 * used direcly as the backing storage for the matrix. Note the use of `STORAGE_USER_OWNED`
 * for the SimpleMatrix's storage policy template parameter.
 *
 * In c++11, a template alias is available for this construction, called `WrapperMatrix`:
 *
 *     WrapperMatrix<double,4,4> wmtx(myRowMajorArray);
 * 
 * For more on matrices, see the @link matrix matrix module documentation@endlink.
 */
template <typename T, index_t M, index_t N, StoragePolicy P>
#ifdef PARSING_DOXYGEN
class SimpleMatrix : public detail::WriteableMatrixBase<T,M,N, SimpleMatrix<T,M,N,P> > {
#else
class FlatMatrixBase : public detail::WriteableMatrixBase<T,M,N, FlatMatrixBase<T,M,N,P> > {
#endif

    typename GenericStorage<T, M*N, P>::type data;
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

#ifdef PARSING_DOXYGEN

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
     * This constructor is not available if the storage policy is `STORAGE_USER_OWNED`.
     * 
     * @param nrows Number of rows in the matrix.
     * @param ncols Number of columns in the matrix.
     */
    explicit SimpleMatrix(index_t nrows, index_t ncols) {}


    /**
     * @brief Construct a matrix of size `(nrows x ncols)`, initialized with
     * `src_data`.
     * 
     * For matrices with dynamic dimension, a size argument is **required** for 
     * that dimension. Size arguments will be ignored for dimensions that are 
     * statically-sized.
     *
     * If the storage policy is `STORAGE_USER_OWNED`, `src_data` will be used directly
     * as the backing storage, and its lifetime must exceed the lifetime of this matrix.
     * 
     * @param src_data Array of nrows * ncols elements, in row-major order.
     * @param nrows Number of rows in the matrix.
     * @param ncols Number of columns in the matrix.
     */
    explicit SimpleMatrix(T* src_data, index_t nrows, index_t ncols) {}

    /**
     * Construct and initialize this matrix with the contents of another.
     *
     * This constructor is not available if the storage policy is `STORAGE_USER_OWNED`.
     * 
     * @tparam Mx A matrix type with agreeing dimension. 
     * @param mtx Matrix containing source elements.
     */
    template <typename Mx> SimpleMatrix(const Mx &mtx) {}

#else

protected:

    explicit FlatMatrixBase(
                 index_t nrows, 
                 index_t ncols,
                 const T* src_data) : 
                     data(nrows * ncols, src_data) {
        Dimension<M>::set(n_rows, nrows);
        Dimension<N>::set(n_cols, ncols);
    }
    
    // const overload of above.
    // needed for user owned storage.
    explicit FlatMatrixBase(
                 index_t nrows, 
                 index_t ncols,
                 T* src_data) : 
                     data(nrows * ncols, src_data) {
        Dimension<M>::set(n_rows, nrows);
        Dimension<N>::set(n_cols, ncols);
    }

    // Undefined behavior if storage backing is user-owned. However,
    // SimpleMatrices with user-owned storage will never try to call this.
    FlatMatrixBase(index_t nrows, index_t ncols) :
            data(nrows * ncols) {
        Dimension<M>::set(n_rows, nrows);
        Dimension<N>::set(n_cols, ncols);
        setIdentity();
    }

public:

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
    typename boost::enable_if<boost::is_scalar<U>, SimpleMatrix<T,M,N,P>&>::type
    operator*=(U k) {
        T *p = data.get();
        for (index_t i = 0; i < rows() * cols(); i++) {
            p[i] *= k;
        }
        return *(static_cast<SimpleMatrix<T,M,N,P>*>(this));
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
        T* elem = data.get() + r * cols() + c;
        *elem = val;
        return *elem;
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

}; // class FlatMatrixBase


#ifndef PARSING_DOXYGEN
} // namespace detail
#endif


///////////////////// SimpleMatrix /////////////////////


// todo: really double check the rows/cols default argument business.


#ifndef PARSING_DOXYGEN


template <typename T, index_t M, index_t N, StoragePolicy P>
class SimpleMatrix : public detail::FlatMatrixBase<T,M,N,P> {

private:

    typedef detail::FlatMatrixBase<T,M,N,P> parent_t;

public:

    // I would not expect this trickery to work, but it does!
    // this mandates row, column arguments for dynamic matrices.
    // arguments are basically ignored if matrix is static.
    explicit SimpleMatrix(
            index_t nrows=detail::DefinedIf<M != DYNAMIC_DIM, M>::value, 
            index_t ncols=detail::DefinedIf<N != DYNAMIC_DIM, N>::value) : 
                parent_t(nrows, ncols) {}

    explicit SimpleMatrix(
            const T* src_data,
            index_t nrows=detail::DefinedIf<M != DYNAMIC_DIM, M>::value, 
            index_t ncols=detail::DefinedIf<N != DYNAMIC_DIM, N>::value) : 
                parent_t(nrows, ncols, src_data) {}

    template <typename Mx>
    SimpleMatrix(const Mx &mtx,
                 typename boost::enable_if_c<
                    (detail::LinalgDimensionMatch<SimpleMatrix<T,M,N>, Mx>::val),
                    void*>::type dummy=0):
                        parent_t(mtx.rows(), mtx.cols()) {
        detail::_mtxcopy(this, mtx);
    }


};  // class SimpleMatrix <...>



///////////////////// SimpleMatrix "wrapper" specialization /////////////////////


template <typename T, index_t M, index_t N>
class SimpleMatrix<T,M,N,STORAGE_USER_OWNED> : public detail::FlatMatrixBase<T,M,N,STORAGE_USER_OWNED> {

private:

    typedef detail::FlatMatrixBase<T,M,N,STORAGE_USER_OWNED> parent_t;

public:

    explicit SimpleMatrix(
            T* src_data,
            index_t nrows=detail::DefinedIf<M != DYNAMIC_DIM, M>::value, 
            index_t ncols=detail::DefinedIf<N != DYNAMIC_DIM, N>::value) : 
                parent_t(nrows, ncols, src_data) {}


};  // class SimpleMatrix <..., USER_OWNER>


#endif


} // end namespace geom

#endif /* SIMPLEMATRIX_H_ */
