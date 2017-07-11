/*
 * MatrixBase.h
 *
 *  Created on: Apr 17, 2013
 *      Author: tbabb
 */

#ifndef NEWMATRIX_H_
#define NEWMATRIX_H_

#include <boost/static_assert.hpp>
#include <algorithm>
#include <geomc/linalg/Vec.h>
#include <geomc/shape/Rect.h>
#include <geomc/linalg/LinalgTypes.h>
#include <geomc/linalg/mtxdetail/MatrixGlue.h>
#include <geomc/linalg/mtxdetail/MatrixIterator.h>

/*
 * Matrix typedefs:
 * 
 *       elem_t           (T)
 *       reference        (T& or AssignmentProxy)
 *       iterator         (T* or SubsetIterator)
 * const_iterator         (const T* or SubsetIterator<const T>)
 *       row_iterator     (T* or SliceIterator or StrideSliceIterator)
 * const_row_iterator     (const T* or SliceIterator<const> or StrideSliceIterator<const>)
 *       col_iterator     (T* or SliceIterator or StrideSliceIterator)
 * const_col_iterator     (const T* or SliceIterator<const> or StrideSliceIterator<const>)
 *       region_iterator  (SubsetIterator)
 * const_region_iterator  (SubsetIterator<const>)
 * 
 * constant iterators cannot be written through. 
 *   - They accomplish this by returning `const T` objects, which are not lvals.
 *   - Non-constant iterators, by contrast, return references (or proxy references) to T objects.
 *   - Constant iterators shall call the const version of m.get(), which by nature
 *     does not return writable references. This is done by passing `const M` as 
 *     the matrix type to their templates.
 */

namespace geom {
namespace detail {

// here we use the curiously recurring template pattern
// to provide "template inheritance", i.e. the ability
// for a whole host of distinct template classes to statically 
// inherit functionality from a "base" template, without virtual
// inheritance and all the runtime costs associated with it.

// this creates a bizarre (but useful) case where the base
// inherits functionality from the derived class.

template <typename T, index_t M, index_t N, typename Derived>
class MatrixBase {

    inline const Derived* derived() const { return static_cast<const Derived*>(this); }
    
    typedef typename detail::_ImplMtxConstRowtype<Derived,T>::const_row_iterator derived_const_row_iterator;
    
public:
    
    /// @brief Row dimension template parameter.
    static const index_t ROWDIM = M;
    /// @brief Column dimension template parameter
    static const index_t COLDIM = N;
    /// @brief Element type
    typedef T elem_t;
    
    /// @brief Read-only row-major iterator over matrix elements.
    typedef MtxSubsetIterator<const Derived,const elem_t> const_iterator;
    /// @brief Read-only row-major iterator over the matrix elements in a rectangular region.
    typedef const_iterator                                const_region_iterator;
    /// @brief Read-only iterator over the elments of a row
    typedef MtxRowIterator<const Derived,const elem_t>    const_row_iterator;
    /// @brief Read-only iterator over the elements of a column
    typedef MtxColIterator<const Derived,const elem_t>    const_col_iterator;
    
    typedef Storage<storage_id_t, _ImplStorageObjCount<Derived>::count> storagebuffer_t;

    // This is accessible from the derived class(es).
    // We need this so we can check boost::is_base_of<>; because
    // the base we are testing has Derived as a template parameter!
    // We need some way of know which class it was that
    // was passed to the curiously-recurring base.
    // I dub this "rectocranial inversion template pattern".
    // Obligatory note that c++ is a stinky pile of poo.
    typedef Derived recurring_t;
    
    /**
     * @return Number of rows in the matrix.
     */
    inline index_t rows() const {
        // "downward" inheritance is weird.
        return derived()->rows();
    }
    
    /**
     * @return Number of columns in the matrix.
     */
    inline index_t cols() const {
        return derived()->cols();
    }
    
    /**
     * @param i Index of row (zero-indexed)
     * @return A const iterator over the elements of row `i`
     */
    inline derived_const_row_iterator operator[](index_t i) const {
        return derived()->row(i);
    }
    
    /**
     * @param i Index of row (zero-indexed)
     * @return A const iterator over the elements of row `i`
     */
    inline const_row_iterator row(index_t i) const {
        return const_row_iterator(derived(), i);
    }
    
    /**
     * @param i Index of column (zero-indexed)
     * @return A const iterator over the elements of column `i`
     */
    inline const_col_iterator col(index_t i) const {
        return const_col_iterator(derived(), i);
    }

    /**
     * Get the element at `(row, col)`.
     * @param row Zero-indexed row coordinate
     * @param col Zero-indexed column coordinate
     * @return The element at `(row, col)`
     */
    inline elem_t get(index_t row, index_t col) const {
#ifdef GEOMC_MTX_CHECK_BOUNDS
            if (row >= rows() || col >= cols() || row < 0 || col < 0) {
                throw std::out_of_range();
            }
#endif
        return derived()->get(row,col);
    }
    
    /**
     * @return A read-only random-access row major-ordered iterator over the elements of this matrix,
     * pointing to the element at (0,0).
     */
    inline const_iterator begin() const {
        return const_iterator(derived());
    }
    
    /**
     * @return A read-only random-access row major-ordered iterator over the elements of this matrix,
     * pointing to the element just beyond the last element in the lower right corner.
     */
    inline const_iterator end() const {
        // TODO: delegate calculation.
        MatrixCoord p = MatrixCoord(rows(), 0);
        return const_iterator(derived(), p);
    }
    
    /**
     * @param r The zero-indexed region to iterate over. The upper extreme coordinates
     * represent the index just beyond the last element to be iterated over.
     * 
     * @return A read-only, random-access, row-major iterator over the elements in region `r`,
     * pointing at the first element in the region (upper left corner).
     */
    inline const_region_iterator region_begin(const MatrixRegion &r) const {
        return const_region_iterator(derived(), r);
    }
    
    /**
     * @param r The zero-indexed region to iterate over. The upper extreme coordinates
     * represent the index just beyond the last element to be iterated over.
     * 
     * @return A read-only, random-access, row-major iterator over the elements in region `r`,
     * pointing at the element just beyond the last element in the region.
     */
    inline const_region_iterator region_end(const MatrixRegion &r) const {
        return const_region_iterator(derived(), r).end();
    }

    // this and its related functions assist in determining
    // whether two matrices alias memory.
    inline storagebuffer_t getStorageIDBuffer() const {
        return storagebuffer_t(derived()->getStorageIDCount());
    }
    
};

/************************************
 * WriteableMatrixBase              *
 ************************************/

template <typename T, index_t M, index_t N, typename Derived>
class WriteableMatrixBase : public MatrixBase<T,M,N,Derived> {
    
    typedef MatrixBase<T,M,N,Derived> parent_t;
    typedef typename detail::_ImplMtxReftype<Derived,T>::reference    derived_reference;
    typedef typename detail::_ImplMtxRowtype<Derived,T>::row_iterator derived_row_iterator;
    
    inline Derived* derived() { return static_cast<Derived*>(this); }
    
public:
    
    /// @brief Reference to element type.
    typedef MtxAssignmentProxy<Derived,T>                reference;
    /// @brief Writeable row-major iterator
    typedef MtxSubsetIterator<Derived,derived_reference> iterator;   
    /// @brief Writeable row-major region iterator
    typedef iterator                                     region_iterator;
    /// @brief Writeable iterator over row elements
    typedef MtxRowIterator<Derived,derived_reference>    row_iterator;
    /// @brief Writeable iterator over column elements
    typedef MtxColIterator<Derived,derived_reference>    col_iterator;

    using parent_t::operator[];
    using parent_t::get;
    using parent_t::row;
    using parent_t::col;
    using parent_t::begin;
    using parent_t::end;
    using parent_t::region_begin;
    using parent_t::region_end;
    
    /**
     * @param i Index of row (zero-indexed)
     * @return A writeable iterator over the elements of row `i`.
     */
    inline derived_row_iterator operator[](index_t i) {
        return derived()->row(i);
    }
    
    /**
     * @param i Index of row (zero-indexed)
     * @return A writeable iterator over the elements of row `i`.
     */
    inline row_iterator row(index_t i) {
        return row_iterator(derived(), i);
    }
    
    /**
     * @param i Index of column (zero-indexed)
     * @return A writeable iterator over the elements of column `i`.
     */
    inline col_iterator col(index_t i) {
        return col_iterator(derived(), i);
    }
    
    /**
     * Get the element at `(row, col)`.
     * @param row Zero-indexed row coordinate
     * @param col Zero-indexed column coordinate
     * @return A reference to the element at `(row, col)`
     */
    inline derived_reference get(index_t row, index_t col) {
#ifdef GEOMC_MTX_CHECK_BOUNDS
            if (row >= derived()->rows() || col >= derived()->cols() || row < 0 || col < 0) {
                throw std::out_of_range();
            }
#endif
        return derived_reference(derived(), row, col);
    }
    
    /**
     * Set the element at `(row, col)` to `val`.
     * @param row Zero-indexed row coordinate
     * @param col Zero-indexed column coordinate
     * @param val New value of element at `(row, col)`
     * @return A reference to the element at `(row, col)`, for convenience.
     */
    inline derived_reference set(index_t row, index_t col, T val) {
#ifdef GEOMC_MTX_CHECK_BOUNDS
            if (row >= derived()->rows() || col >= derived()->cols() || row < 0 || col < 0) {
                throw std::out_of_range();
            }
#endif
        return derived()->set(row,col, val);
    }
    
    /**
     * @return A writeable, random-access, row-major iterator over the elements of this matrix,
     * pointing to the element at (0,0).
     */
    inline iterator begin() {
        return iterator(derived());
    }
    
    /**
     * @return A writeable, random-access, row-major iterator over the elements of this matrix,
     * pointing to the element just beyond the last element in the lower right corner.
     */
    inline iterator end() {
        // TODO: delegate calculation.
        MatrixCoord p = MatrixCoord(derived()->rows(), 0);
        return iterator(derived(), p);
    }
    
    /**
     * @param r The zero-indexed region to iterate over. The upper extreme coordinates
     * represent the index just beyond the last element to be iterated over.
     * 
     * @return A writeable, random-access, row-major iterator over the elements in region `r`,
     * pointing at the first element in the region (upper left corner).
     */
    inline region_iterator region_begin(const MatrixRegion &r) {
        return region_iterator(derived(), r);
    }
    
    /**
     * @param r The zero-indexed region to iterate over. The upper extreme coordinates
     * represent the index just beyond the last element to be iterated over.
     * 
     * @return A writeable, random-access, row-major iterator over the elements in region `r`,
     * pointing at the element just beyond the last element in the region.
     */
    inline region_iterator region_end(const MatrixRegion &r) {
        return region_iterator(derived(), r).end();
    }
    
}; // class WriteableMatrixBase

}; // namespace detail
}; // namespace geom


#endif /* NEWMATRIX_H_ */
