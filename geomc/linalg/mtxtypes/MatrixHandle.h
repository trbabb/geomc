/*
 * MatrixHandle.h
 *
 *  Created on: Apr 28, 2013
 *      Author: tbabb
 */

#ifndef GENERALMATRIX_H_
#define GENERALMATRIX_H_

#include <geomc/linalg/mtxdetail/MatrixBase.h>

namespace geom {

namespace detail {


    // const T because m[r][c].mutate() may not make sense. m[r][c] = X is allowed, however.
    template <typename T>
    struct _ImplMtxReftype<MatrixHandle<T>, T> {
        typedef detail::MtxAssignmentProxy<MatrixHandle<T>, const T> reference;
    };
    
    // storage count is dynamic
    template <typename T>
    struct _ImplStorageObjCount<MatrixHandle<T> > {
        const static index_t count;
    };
    
    template <typename T>
    const index_t _ImplStorageObjCount<MatrixHandle<T> >::count = DYNAMIC_DIM;

    
}; // end namespace detail

/** @ingroup matrix 
 *  @brief A generic matrix class which can hold references to all other
 *  matrix types.
 */
template <typename T>
class MatrixHandle : public detail::WriteableMatrixBase<T,DYNAMIC_DIM, DYNAMIC_DIM, MatrixHandle<T> > {
public:
    typedef typename detail::_ImplMtxReftype<MatrixHandle<T>,T>::reference reference;
    
    using detail::WriteableMatrixBase<T,DYNAMIC_DIM,DYNAMIC_DIM, MatrixHandle<T> >::get;
    
    MatrixHandle() {}
    virtual ~MatrixHandle() {}
    
    //////////// methods ////////////
    
    virtual index_t   rows() const = 0;
    virtual index_t   cols() const = 0;
    virtual T         operator()(index_t r, index_t c) const = 0;
    virtual reference set(index_t r, index_t c, T val) = 0;
    virtual void      set_identity() = 0;
    virtual void      set_zero() = 0;
    
    virtual void      get_storage_tokens(storage_token_t* buf) const = 0;
    virtual index_t   storage_token_count() const = 0;
    
};

template <typename M>
class MatrixWrapper : public MatrixHandle<typename M::elem_t> {
public:
    
    typedef typename M::elem_t elem_t;
    typedef typename MatrixHandle<elem_t>::reference reference;
    
    using MatrixHandle<elem_t>::get;
    
    M *m;
    MatrixWrapper(M *m):m(m) {}
    virtual ~MatrixWrapper() {}
    
    //////////// methods ////////////
    
    virtual elem_t operator()(index_t r, index_t c) const {
        return (*m)(r, c);
    }
    
    virtual reference set(index_t r, index_t c, elem_t val) {
        m->set(r, c, val);
        return reference(this,r,c);
    }
    
    virtual void set_identity() {
        m->set_identity();
    }
    
    virtual void set_zero() {
        m->set_zero();
    }
    
    virtual index_t rows() const {
        return m->rows();
    }
    
    virtual index_t cols() const {
        return m->cols();
    }
    
    virtual void get_storage_tokens(storage_token_t* buf) const {
        m->get_storage_tokens(buf);
    }
    
    virtual index_t storage_token_count() const {
        return m->storage_token_count();
    }
};

}; // end namespace geom


#endif /* GENERALMATRIX_H_ */
