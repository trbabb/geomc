/*
 * GeneralMatrix.h
 *
 *  Created on: Apr 28, 2013
 *      Author: tbabb
 */

#ifndef GENERALMATRIX_H_
#define GENERALMATRIX_H_

#include "linalg/mtxdetail/MatrixBase.h"

namespace geom {

namespace detail {


    // const T because m[r][c].mutate() may not make sense. m[r][c] = X is allowed, however.
    template <typename T>
    struct _ImplMtxReftype<Matrix<T>, T> {
        typedef detail::MtxAssignmentProxy<Matrix<T>, const T> reference;
    };
    
    // storage count is dynamic
    template <typename T>
    struct _ImplStorageObjCount<Matrix<T> > {
        const static index_t count;
    };
    
    template <typename T>
    const index_t _ImplStorageObjCount<Matrix<T> >::count = DYNAMIC_DIM;

    
}; // end namespace detail


template <typename T>
class Matrix : public detail::WriteableMatrixBase<T,DYNAMIC_DIM, DYNAMIC_DIM, Matrix<T> > {
public:
    typedef typename detail::_ImplMtxReftype<Matrix<T>,T>::reference reference;
    
    using detail::WriteableMatrixBase<T,DYNAMIC_DIM,DYNAMIC_DIM, Matrix<T> >::get;
    
    Matrix() {}
    virtual ~Matrix() {}
    
    //////////// methods ////////////
    
    virtual index_t   rows() const = 0;
    virtual index_t   cols() const = 0;
    virtual T         get(index_t r, index_t c) const = 0;
    virtual reference set(index_t r, index_t c, T val) = 0;
    virtual void      setIdentity() = 0;
    virtual void      setZero() = 0;
    
    virtual void      getStorageIDs(storage_id_t *buf) const = 0;
    virtual index_t   getStorageIDCount() const = 0;
    
};

template <typename M>
class MatrixHandle : public Matrix<typename M::elem_t> {
public:
    
    typedef typename M::elem_t elem_t;
    typedef typename Matrix<elem_t>::reference reference;
    
    using Matrix<elem_t>::get;
    
    M *m;
    MatrixHandle(M *m):m(m) {}
    virtual ~MatrixHandle() {}
    
    //////////// methods ////////////
    
    virtual elem_t get(index_t r, index_t c) const {
        return m->get(r, c);
    }
    
    virtual reference set(index_t r, index_t c, elem_t val) {
        m->set(r, c, val);
        return reference(this,r,c);
    }
    
    virtual void setIdentity() {
        m->setIdentity();
    }
    
    virtual void setZero() {
        m->setZero();
    }
    
    virtual index_t rows() const {
        return m->rows();
    }
    
    virtual index_t cols() const {
        return m->cols();
    }
    
    virtual void getStorageIDs(storage_id_t *buf) const {
        m->getStorageIDs(buf);
    }
    
    virtual index_t getStorageIDCount() const {
        return m->getStorageIDCount();
    }
};

}; // end namespace geom


#endif /* GENERALMATRIX_H_ */
