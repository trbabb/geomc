/*
 * NewMatrix.h
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

#if GEOMC_MTX_CHECK_BOUNDS
    #include <geomc/GeomException.h>
#endif

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

// TODO: conversion operator for simplematrix for all other matrix types
//       so that `SimpleMtx m = mtx * mtx` will always work.

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

    inline const Derived *derived() const { return static_cast<const Derived*>(this); }
    
    typedef typename detail::_ImplMtxConstRowtype<Derived,T>::const_row_iterator derived_const_row_iterator;
    
public:
    
    static const index_t ROWDIM = M;
    static const index_t COLDIM = N;
    typedef T elem_t;
    
    typedef MtxSubsetIterator<const Derived,const elem_t> const_iterator;
    typedef const_iterator                                const_region_iterator;
    typedef MtxRowIterator<const Derived,const elem_t>    const_row_iterator;
    typedef MtxColIterator<const Derived,const elem_t>    const_col_iterator;
    
    typedef Storage<storage_id_t, _ImplStorageObjCount<Derived>::count> storagebuffer_t;
    
    inline index_t rows() const {
        // "downward" inheritance is weird.
        return derived()->rows();
    }
    
    inline index_t cols() const {
        return derived()->cols();
    }
    
    inline derived_const_row_iterator operator[](index_t i) const {
        return derived()->row(i);
    }
    
    inline const_row_iterator row(index_t i) const {
        return const_row_iterator(derived(), i);
    }
    
    inline const_col_iterator col(index_t i) const {
        return const_col_iterator(derived(), i);
    }

    inline elem_t get(index_t row, index_t col) const {
        #ifdef GEOMC_MTX_CHECK_BOUNDS
            if (row >= rows() || col >= cols() || row < 0 || col < 0) {
                throw std::out_of_range();
            }
        #endif
        return derived()->get(row,col);
    }
    
    inline const_iterator begin() const {
        return const_iterator(derived());
    }
    
    inline const_iterator end() const {
        // TODO: delegate calculation.
        MatrixCoord p = MatrixCoord(rows(), 0);
        return const_iterator(derived(), p);
    }
    
    inline const_region_iterator region_begin(const MatrixRegion &r) const {
        return const_region_iterator(derived(), r);
    }
    
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
    
    typedef MtxAssignmentProxy<Derived,T>                reference;
    typedef MtxSubsetIterator<Derived,derived_reference> iterator;   
    typedef iterator                                     region_iterator;
    typedef MtxRowIterator<Derived,derived_reference>    row_iterator;
    typedef MtxColIterator<Derived,derived_reference>    col_iterator;

    using parent_t::operator[];
    using parent_t::get;
    using parent_t::row;
    using parent_t::col;
    using parent_t::begin;
    using parent_t::end;
    using parent_t::region_begin;
    using parent_t::region_end;
    
    inline derived_row_iterator operator[](index_t i) {
        return derived()->row(i);
    }
    
    inline row_iterator row(index_t i) {
        return row_iterator(derived(), i);
    }
    
    inline col_iterator col(index_t i) {
        return col_iterator(derived(), i);
    }

    inline derived_reference get(index_t row, index_t col) {
        #ifdef GEOMC_MTX_CHECK_BOUNDS
            if (row >= derived()->rows() || col >= derived()->cols() || row < 0 || col < 0) {
                throw std::out_of_range();
            }
        #endif
        return derived_reference(derived(), row, col);
    }
    
    inline derived_reference set(index_t row, index_t col) {
        #ifdef GEOMC_MTX_CHECK_BOUNDS
            if (row >= derived()->rows() || col >= derived()->cols() || row < 0 || col < 0) {
                throw std::out_of_range();
            }
        #endif
        return derived()->set(row,col);
    }
    
    inline iterator begin() {
        return iterator(derived());
    }
    
    inline iterator end() {
        // TODO: delegate calculation.
        MatrixCoord p = MatrixCoord(derived()->rows(), 0);
        return iterator(derived(), p);
    }
    
    inline region_iterator region_begin(const MatrixRegion &r) {
        return region_iterator(derived(), r);
    }
    
    inline region_iterator region_end(const MatrixRegion &r) {
        return region_iterator(derived(), r).end();
    }
    
}; // class WriteableMatrixBase

}; // namespace detail
}; // namespace geom


#endif /* NEWMATRIX_H_ */
