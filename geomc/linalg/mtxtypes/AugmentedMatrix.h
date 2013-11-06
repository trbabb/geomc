/*
 * AugmentedMatrix.h
 *
 *  Created on: Apr 28, 2013
 *      Author: tbabb
 */

#ifndef AUGMENTEDMATRIX_H_
#define AUGMENTEDMATRIX_H_

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <geomc/linalg/mtxdetail/MatrixBase.h>

namespace geom {

namespace detail {

    // const T because m[r][c].mutate() does not make sense. m[r][c] = X is allowed, however.
    template <typename Ma, typename Mb>
    struct _ImplMtxReftype<AugmentedMatrix<Ma,Mb>, typename Ma::elem_t> {
        typedef detail::MtxAssignmentProxy<AugmentedMatrix<Ma,Mb>, const typename Ma::elem_t> reference;
    };
    
    // storage count
    template <typename Ma, typename Mb>
    struct _ImplStorageObjCount<AugmentedMatrix<Ma,Mb> > {
        const static index_t count;
    };
    
    // if either matrix has dynamic count, then we must have dynamic count
    template <typename Ma, typename Mb>
    const index_t _ImplStorageObjCount<AugmentedMatrix<Ma,Mb> >::count = 
            (_ImplStorageObjCount<Ma>::count == DYNAMIC_DIM || _ImplStorageObjCount<Mb>::count == DYNAMIC_DIM) ?
                    (DYNAMIC_DIM) : 
                    (_ImplStorageObjCount<Ma>::count + _ImplStorageObjCount<Mb>::count);
    
}; // end namespace detail


// TODO: The reference type mutability may not match
// TODO: the mutability of the child classes may not match
// TODO: ownership of child ptrs?

/** @ingroup matrix 
 *  @brief A matrix which wraps two side-by-side sub-matrices.
 */
template <typename Ma, typename Mb>
class AugmentedMatrix : public detail::WriteableMatrixBase<typename Ma::elem_t, DYNAMIC_DIM, DYNAMIC_DIM, AugmentedMatrix<Ma,Mb> > {
    
    Ma* mtx_a;
    Mb* mtx_b;
    index_t n_rows;
    index_t n_cols;
    
    //TODO: figure out why no worky
    //BOOST_STATIC_ASSERT((boost::is_same<typename Ma::elem_t, typename Mb::elem_t>::value));
    
public:
    
    typedef typename Ma::elem_t elem_t;
    typedef typename detail::_ImplMtxReftype<AugmentedMatrix<Ma,Mb>,elem_t>::reference reference;
    
    using detail::WriteableMatrixBase<typename Ma::elem_t, DYNAMIC_DIM, DYNAMIC_DIM, AugmentedMatrix<Ma,Mb> >::get;
    
    AugmentedMatrix(Ma *m_a, Mb *m_b):
        mtx_a(m_a), 
        mtx_b(m_b),
        n_rows(std::min(mtx_a->rows(), mtx_b->rows())),
        n_cols(mtx_a->cols() + mtx_b->cols()) {
        // TODO: bounds check.
    }
    
    inline index_t rows() const {
        return n_rows;
    }
    
    inline index_t cols() const {
        return n_cols;
    }
    
    inline elem_t get(index_t r, index_t c) const {
#ifdef GEOMC_MTX_CHECK_BOUNDS
            if (r < 0 or r >= rows() or c < 0 or c >= cols()) {
                throw std::out_of_range("matrix coordinates");
            }
#endif
        if (c < mtx_a->cols()) {
            return ((const Ma*)mtx_a)->get(r, c);
        } else {
            return ((const Mb*)mtx_b)->get(r, c - mtx_a->cols());
        }
    }
    
    inline reference set(index_t r, index_t c, elem_t val) {
        if (c < mtx_a->cols()) {
            mtx_a->set(r, c, val);
        } else {
            mtx_b->set(r, c - mtx_a->cols(), val);
        }
        return reference(this, r, c);
    }

    void setIdentity() {
        setZero();
        for (index_t i = 0; i < std::min(rows(), cols()); i++) {
            set(i,i,1);
        }
    }
    
    inline void setZero() {
        mtx_a->setZero();
        mtx_b->setZero();
    }
    
    inline void getStorageIDs(storage_id_t *buf) const {
        mtx_a->getStorageIDs(buf);
        mtx_b->getStorageIDs(buf + mtx_a->getStorageIDCount());
    }
    
    inline index_t getStorageIDCount() const {
        return mtx_a->getStorageIDCount() + mtx_b->getStorageIDCount();
    }
};

}; // end namespace geom


#endif /* AUGMENTEDMATRIX_H_ */
