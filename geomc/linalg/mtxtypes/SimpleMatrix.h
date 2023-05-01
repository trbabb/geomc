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

#include <type_traits>
#include <geomc/Hash.h>
#include <geomc/Storage.h>
#include <geomc/linalg/mtxdetail/MatrixLayout.h>
#include <geomc/linalg/mtxdetail/MatrixBase.h>
#include <geomc/linalg/mtxdetail/MatrixCopy.h>


/* 
    what would we need for templated layout?

    - tell us the row iterators
    - tell us the column iterators
    - provide get() --> {const,} iterator
    - tell us how create storage
    - tell us how to fill(?)
    
    Would need to update generic row and column iterators to
    iterate over the whole matrix, not just their own column.
    This would mean that row() or col() would behave the same
    way regardless of layout.
    
    This would allow us to more easily treat vectors as columns
    directly; e.g. for linear solve, tracing triangle, orthogonalizaing...
    
    Might be even better to abstract the memory into a Grid<...>.
*/

namespace geom {


///////////////////// Details /////////////////////


namespace detail {

    // allow myself to introduce myself:
    // I am a curiously recurring class, and I have to declare myself to declare myself.
    template <typename T, index_t M, index_t N, MatrixLayout Lyt, StoragePolicy P>
    class FlatMatrixBase;

    // specializations of the iterator type helper classes
    // tell the base class which iterators to use:

    template <typename T, index_t M, index_t N, MatrixLayout Lyt, StoragePolicy P>
    struct _ImplMtxReftype<FlatMatrixBase<T,M,N,Lyt,P>,T> {
        typedef T& reference;
    };
    
    template <typename T, index_t M, index_t N, MatrixLayout Lyt, StoragePolicy P>
    struct _ImplMtxRowtype<FlatMatrixBase<T,M,N,Lyt,P>,T> {
        typedef typename FlatMatrixLayout<T,Lyt>::row_iterator row_iterator;
    };
    
    template <typename T, index_t M, index_t N, MatrixLayout Lyt, StoragePolicy P>
    struct _ImplMtxConstRowtype<FlatMatrixBase<T,M,N,Lyt,P>,T> {
        typedef typename FlatMatrixLayout<T,Lyt>::const_row_iterator const_row_iterator;
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
 * @tparam Lyt Layout of the underlying array (ROW_MAJOR or COL_MAJOR). ROW_MAJOR is the default.
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
 * <li>If the StoragePolicy is `STORAGE_WRAPPED`, then the matrix will wrap a user-provided
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
 *     SimpleMatrix<double,4,4,STORAGE_WRAPPED> mtx(myRowMajorArray);
 *
 * In the above example, no array duplication will occur, and `myRowMajorArray` is 
 * used direcly as the backing storage for the matrix. Note the use of `STORAGE_WRAPPED`
 * for the SimpleMatrix's storage policy template parameter.
 *
 * In c++11, a template alias is available for this construction, called `WrapperMatrix`:
 *
 *     WrapperMatrix<double,4,4> wmtx(myRowMajorArray);
 * 
 * For more on matrices, see the @link matrix matrix module documentation@endlink.
 */
template <typename T, index_t M, index_t N, MatrixLayout Lyt, StoragePolicy P>
#ifdef PARSING_DOXYGEN
class SimpleMatrix : public detail::WriteableMatrixBase<T,M,N, SimpleMatrix<T,M,N,Lyt,P> > {
#else
class FlatMatrixBase : public detail::WriteableMatrixBase<T,M,N, FlatMatrixBase<T,M,N,Lyt,P> > {
#endif

    typename GenericStorage<T, M * N, P>::type data;
    typename Dimension<M>::storage_t n_rows; // this contortion prevents allocating storage 
    typename Dimension<N>::storage_t n_cols; // for a dimension size when it is not dynamic.
    
    static_assert(M >= 0, "Row dimension must not be negative");
    static_assert(N >= 0, "Column dimension must not be negative");
    
public:
    
    // MatrixBase types
    /**
     * @brief Read-only row-major iterator over the elements of the matrix.
     * This is a `T*` if the matrix is row-major; a proxy otherwise.
     */ 
    typedef typename FlatMatrixLayout<T,Lyt>::const_row_iterator const_iterator;
    /**
     * @brief Read-only iterator over the elements of a row.
     * This is a `T*` if the matrix is row-major; a proxy otherwise.
     */
    typedef typename FlatMatrixLayout<T,Lyt>::const_row_iterator const_row_iterator;
    /**
     * @brief Read-only iterator over the elements of a column.
     * This is a `T*` if the matrix is column-major; a proxy otherwise.
     */
    typedef typename FlatMatrixLayout<T,Lyt>::const_col_iterator const_col_iterator;
    
    // WriteableMatrixBase types
    /// @brief Reference to element type.
    typedef T& reference;
    /**
     * @brief Writeable row-major iterator over the elements of the matrix.
     * This is a `T*` if the matrix is row-major; a proxy otherwise.
     */ 
    typedef typename FlatMatrixLayout<T,Lyt>::row_iterator iterator;
    /**
     * @brief Writeable iterator over the elements of a row.
     * This is a `T*` if the matrix is row-major; a proxy otherwise.
     */
    typedef typename FlatMatrixLayout<T,Lyt>::row_iterator row_iterator;
    /**
     * @brief Writeable iterator over the elements of a column.
     * This is a `T*` if the matrix is column-major; a proxy otherwise.
     */
    typedef typename FlatMatrixLayout<T,Lyt>::col_iterator col_iterator;
    
    /// The data layout of this matrix.
    static constexpr MatrixLayout Layout = Lyt;

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
     * This constructor is not available if the storage policy is `STORAGE_WRAPPED`.
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
     * If the storage policy is `STORAGE_WRAPPED`, `src_data` will be used directly
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
     * This constructor is not available if the storage policy is `STORAGE_WRAPPED`.
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
                 const T* src_data): 
                     data(nrows * ncols, src_data) {
        Dimension<M>::set(n_rows, nrows);
        Dimension<N>::set(n_cols, ncols);
    }
    
    // const overload of above.
    // needed for user owned storage.
    explicit FlatMatrixBase(
                 index_t nrows, 
                 index_t ncols,
                 T* src_data): 
                     data(nrows * ncols, src_data) {
        Dimension<M>::set(n_rows, nrows);
        Dimension<N>::set(n_cols, ncols);
    }

    // Undefined behavior if storage backing is user-owned. However,
    // SimpleMatrices with user-owned storage will never try to call this.
    FlatMatrixBase(index_t nrows, index_t ncols):
            data(nrows * ncols) {
        Dimension<M>::set(n_rows, nrows);
        Dimension<N>::set(n_cols, ncols);
        set_identity();
    }
    
    // copy from another matrix without setting identity first.
    template <typename Mx>
    FlatMatrixBase(const Mx& m):
            data(m.rows() * m.cols()) {
        Dimension<M>::set(n_rows, m.rows());
        Dimension<N>::set(n_cols, m.cols());
        detail::_mtxcopy(this, m);
    }

public:

#endif
    
    /**************************
     * Operators              *
     **************************/
    
    /**
     * Scalar in-place multiplication.
     * @tparam U  A type satisfying `std::is_scalar<>`.
     * @param  k  Scalar value.
     * @return    A reference to `this`. 
     */
#ifdef PARSING_DOXYGEN
    template <typename U> SimpleMatrix<T,M,N>& operator*=(U k) {}
#else
    template <typename U>
    typename std::enable_if<
            std::is_scalar<U>::value,
            SimpleMatrix<T,M,N,Lyt,P>&
        >::type
    operator*=(U k) {
        T *p = data.get();
        for (index_t i = 0; i < rows() * cols(); i++) {
            p[i] *= k;
        }
        return *(static_cast<SimpleMatrix<T,M,N,Lyt,P>*>(this));
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
        return FlatMatrixLayout<T,Lyt>::template row<true>(data.get(), i, rows(), cols());
    }
    
    /**
     * @param i Index of row (zero-indexed)
     * @return A writeable iterator over the elements of row `i`.
     */
    inline row_iterator row(index_t i) {
        return FlatMatrixLayout<T,Lyt>::template row<false>(data.get(), i, rows(), cols());
    }
    
    /**
     * @param i Index of column (zero-indexed)
     * @return A read-only iterator over the elements of row `i`.
     */
    inline const_col_iterator col(index_t i) const {
        return FlatMatrixLayout<T,Lyt>::template col<true>(data.get(), i, rows(), cols());
    }
    
    /**
     * @param i Index of column (zero-indexed)
     * @return A writeable iterator over the elements of row `i`.
     */
    inline col_iterator col(index_t i) {
        return FlatMatrixLayout<T,Lyt>::template col<false>(data.get(), i, rows(), cols());
    }
    
    /**
     * @return A read-only, random-access, row-major iterator over the elements of this matrix,
     * pointing to the element at (0,0).
     */
    inline const_iterator begin() const {
        return FlatMatrixLayout<T,Lyt>::template row<true>(data.get(), 0, rows(), cols());
    }
    
    /**
     * @return A read-only, random-access, row-major iterator over the elements of this matrix,
     * pointing to the just beyond the last element in the matrix.
     */
    inline const_iterator end() const {
        return begin() + (rows() * cols());
    }
    
    /**
     * @return A writeable, random-access, row-major iterator over the elements of this matrix,
     * pointing to the element at (0,0).
     */
    inline iterator begin() {
        return FlatMatrixLayout<T,Lyt>::template row<false>(data.get(), 0, rows(), cols());
    }
    
    /**
     * @return A writeable, random-access, row-major iterator over the elements of this matrix,
     * pointing to the element just beyond the last element in the matrix.
     */
    inline iterator end() {
        return begin() + (rows() * cols());
    }
    
    /**
     * @return A pointer to the bare data array, in the layout policy of this matrix.
     */
    inline T* data_begin() {
        return data.get();
    }
    
    /**
     * @return A pointer to the end of the bare data array.
     */
    inline T* data_end() {
        return data.get() + (rows() * cols());
    }
    
    /**
     * @return A pointer to the read-only bare data array, in the layout policy of this matrix.
     */
    inline const T* data_begin() const {
        return data.get();
    }
    
    /**
     * @return A pointer to the end of the read-only bare data array.
     */
    inline const T* data_end() const {
        return data.get() + (rows() * cols());
    }
    
    
    /**
     * Get the element at `(row, col)`.
     * @param r Zero-indexed row coordinate
     * @param c Zero-indexed column coordinate
     * @return The element at `(row, col)`
     */
    inline T operator()(index_t r, index_t c) const {
        const index_t i = FlatMatrixLayout<T,Lyt>::index(r, c, rows(), cols());
        return data.get()[i];
    }
    
    /**
     * Get the element at `(row, col)`.
     * @param r Zero-indexed row coordinate
     * @param c Zero-indexed column coordinate
     * @return A reference to the element at `(row, col)`
     */
    inline reference operator()(index_t r, index_t c) {
        const index_t i = FlatMatrixLayout<T,Lyt>::index(r, c, rows(), cols());
        return data.get()[i];
    }
    
    /**
     * Set the element at `(row, col)` to `val`.
     * @param r Zero-indexed row coordinate
     * @param c Zero-indexed column coordinate
     * @param val New value of element at `(row, col)`
     * @return A reference to the element at `(row, col)`, for convenience.
     */
    inline reference set(index_t r, index_t c, T val) {
        const index_t i = FlatMatrixLayout<T,Lyt>::index(r, c, rows(), cols());
        T* elem = data.get() + i;
        *elem = val;
        return *elem;
    }
    
    /**
     * Reinterpret this matrix as a transposed version of itself.
     *
     * The returned matrix is backed by this matrix's underlying memory. As such, its 
     * lifetime must not be longer than that of this matrix. Changes to the elements of this
     * matrix will also be reflected in the resultant matrix, and vice versa.
     *
     * This operation is O(1), and fuctions by constructing a new `WrapperMatrix` with opposite
     * layout (i.e., column-major if this matrix is row-major, or row-major if this matrix is
     * column-major).
     */
#ifdef PARSING_DOXYGEN
    WrapperMatrix<T,N,M, (Lyt == ROW_MAJOR) ? COL_MAJOR : ROW_MAJOR> 
    as_transpose() const {
#else
    SimpleMatrix<T,N,M, (Lyt == ROW_MAJOR) ? COL_MAJOR : ROW_MAJOR, STORAGE_WRAPPED> 
    as_transpose() const {
#endif
        return 
        SimpleMatrix<T,N,M, (Lyt == ROW_MAJOR) ? COL_MAJOR : ROW_MAJOR, STORAGE_WRAPPED>(
            data.get(), cols(), rows()
        );
    }
    
    /**
     * Clear this matrix and set its elements to the identity matrix.
     */
    void set_identity() {
        set_zero();
        T *last = data_end();
        for (T *p = data_begin(); p < last; p += cols() + 1) {
            *p = 1;
        }
    }
    
    /**
     * Set all the elements of this matrix to zero.
     */
    inline void set_zero() {
        std::fill(data_begin(), data_end(), 0);
    }
    
    inline void get_storage_tokens(storage_token_t* buf) const {
        const T* p = data.get();
        *buf = {p, p + rows() * cols()};
    }
    
    constexpr index_t storage_token_count() const {
        return 1;
    }

}; // class FlatMatrixBase


#ifndef PARSING_DOXYGEN
} // namespace detail
#endif


///////////////////// SimpleMatrix /////////////////////


// todo: really double check the rows/cols default argument business.


#ifndef PARSING_DOXYGEN


template <typename T, index_t M, index_t N, MatrixLayout Lyt, StoragePolicy P>
class SimpleMatrix : public detail::FlatMatrixBase<T,M,N,Lyt,P> {

private:

    typedef detail::FlatMatrixBase<T,M,N,Lyt,P> parent_t;

public:

    // I would not expect this trickery to work, but it does!
    // this mandates row, column arguments for dynamic matrices.
    // arguments are basically ignored if matrix is static.
    explicit SimpleMatrix(
        // xxx test this
            index_t nrows=detail::DefinedIf<M != DYNAMIC_DIM, M>::value, 
            index_t ncols=detail::DefinedIf<N != DYNAMIC_DIM, N>::value) : 
                parent_t(nrows, ncols) {}

    explicit SimpleMatrix(
            const T* src_data,
            index_t nrows=detail::DefinedIf<M != DYNAMIC_DIM, M>::value, 
            index_t ncols=detail::DefinedIf<N != DYNAMIC_DIM, N>::value) : 
                parent_t(nrows, ncols, src_data) {}
    
    SimpleMatrix(const std::initializer_list<T>& elems):
            parent_t(M,N,elems.begin())
    {
        static_assert(M * N != 0, "initializer list only allowed with static dimensions");
    }

    template <typename Mx>
    SimpleMatrix(
        const Mx &mtx,
        typename std::enable_if
        <
            detail::LinalgDimensionMatch<SimpleMatrix<T,M,N>, Mx>::val,
            void*
        >::type dummy=0):
            parent_t(mtx) {}


};  // class SimpleMatrix <...>



///////////////////// SimpleMatrix "wrapper" specialization /////////////////////


template <typename T, index_t M, index_t N, MatrixLayout Lyt>
class SimpleMatrix<T,M,N,Lyt,STORAGE_WRAPPED>: 
    public detail::FlatMatrixBase<T,M,N,Lyt,STORAGE_WRAPPED>
{

private:

    typedef detail::FlatMatrixBase<T,M,N,Lyt,STORAGE_WRAPPED> parent_t;

public:

    explicit SimpleMatrix(T* src_data):
            parent_t(M, N, src_data) {
        static_assert(M * N != 0,
            "Matrices with dynamic dimensions must specify a runtime size.");
    }
    
    SimpleMatrix(T* src_data, index_t n):
            parent_t(
                M == 0 ? n : M, 
                N == 0 ? n : N, 
                src_data)
    {
        static_assert(std::max(M, N) > 0,
            "Matrix has two dynamic dimensions, but only one runtime size was specified.");
    }
    
    SimpleMatrix(T* src_data, index_t nrows, index_t ncols):
            parent_t(nrows, ncols, src_data) {
        static_assert(N == 0 and M == 0,
            "Matrix has static dimensions, but a runtime dimensions were specified.");
    }


};  // class SimpleMatrix <..., STORAGE_WRAPPED>


#endif


} // end namespace geom


namespace std {

template <typename T, index_t M, index_t N>
struct hash<geom::SimpleMatrix<T,M,N>> {
    size_t operator()(const geom::SimpleMatrix<T,M,N> &m) const {
        return geom::hash_bytes(m.data_begin(), m.rows() * m.cols() * sizeof(T));
    }
};

} // end namespace std

#endif /* SIMPLEMATRIX_H_ */
