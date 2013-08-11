/*
 * SparseMatrix.h
 *
 *  Created on: Apr 28, 2013
 *      Author: tbabb
 */

#ifndef SPARSEMATRIX_H_
#define SPARSEMATRIX_H_

#include <boost/config.hpp>

#ifdef BOOST_NO_CXX11_HDR_UNORDERED_MAP 
    #include <tr1/unordered_map>
#else
    #include <unordered_map>
#endif

#include <geomc/linalg/mtxdetail/MatrixBase.h>

namespace geom {

namespace detail {

    // const T because m[r][c].mutate() may not make sense. m[r][c] = X is allowed, however.
    template <typename T>
    struct _ImplMtxReftype<SparseMatrix<T>,T> {
        typedef detail::MtxAssignmentProxy<SparseMatrix<T>, const T> reference;
    };
    
    // fwd decl
    template <typename T>
    inline void _mtxcopy(geom::SparseMatrix<T> *into, const geom::SparseMatrix<T> &src);
    
}; // end namespace detail

/** @ingroup matrix 
 *  @brief A (still weakly-supported) class which only stores nonzero matrix elements.
 * 
 * Matrix multiplications, and many other operations, are still not yet optimized
 * to take advantage of the sparse represenation. Future versions of the library
 * will likely address this.
 */
template <typename T>
class SparseMatrix : public detail::WriteableMatrixBase<T, DYNAMIC_DIM, DYNAMIC_DIM, SparseMatrix<T> > {
    
#ifdef BOOST_NO_CXX11_HDR_UNORDERED_MAP
    typedef std::tr1::unordered_map<MatrixCoord,T> map_t;
#else
    // use c++11 unordered_map
    typedef std::unordered_map<MatrixCoord,T> map_t;
#endif 
    
    index_t n_rows;
    index_t n_cols;
    map_t elements;
    
public:
    
    typedef typename detail::_ImplMtxReftype<SparseMatrix<T>,T>::reference reference;
    typedef typename map_t::iterator nonzero_iterator;
    
    SparseMatrix(index_t nrows, index_t ncols) : n_rows(nrows), n_cols(ncols) {}
    
    inline index_t rows() const {
        return n_rows;
    }
    
    inline index_t cols() const {
        return n_cols;
    }
    
    // don't hide this; we're overloading based on const-ness.
    using detail::WriteableMatrixBase<T,DYNAMIC_DIM, DYNAMIC_DIM, SparseMatrix<T> >::get;
    
    inline T get(index_t r, index_t c) const {
        #ifdef GEOMC_MTX_CHECK_BOUNDS
            if (r < 0 || r >= rows() || c < 0 || c >= cols()) {
                throw std::out_of_range();
            }
        #endif
        typename map_t::const_iterator i = elements.find(MatrixCoord(r,c));
        if (i == elements.end()) {
            return 0;
        } else {
            return i->second;
        }
    }
    
    inline reference set(index_t r, index_t c, T val) {
        #ifdef GEOMC_MTX_CHECK_BOUNDS
            if (row >= nrows or col >= ncols) {
                throw std::out_of_range("matrix coordinates");
            }
        #endif
        MatrixCoord i = MatrixCoord(r, c);
        if (val == 0) {
            elements.erase(i);
        } else {
            elements.insert(typename map_t::value_type(i, val));
        }
        return reference(this, r, c);
    }
    
    inline nonzero_iterator nonzero_begin() {
        return elements.begin();
    }
    
    inline nonzero_iterator nonzero_end() {
        return elements.end();
    }

    void setIdentity() {
        setZero();
        index_t end = std::min(rows(), cols());
        for (int i = 0; i < end; i++) {
            set(i,i,1);
        }
    }
    
    inline void setZero() {
        elements.clear();
    }
    
    inline void getStorageIDs(storage_id_t *buf) const {
        *buf = &elements;
    }
    
    inline index_t getStorageIDCount() const {
        return 1;
    }

    /******** friends ********/

    friend void detail::_mtxcopy<T>(SparseMatrix<T> *into, const SparseMatrix<T> &src);
};

}; // end namespace geom

#endif /* SPARSEMATRIX_H_ */
