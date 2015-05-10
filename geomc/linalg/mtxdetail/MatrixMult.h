/*
 * MatrixMult.h
 *
 *  Created on: Jun 23, 2013
 *      Author: tbabb
 */

#ifndef MATRIXMULT_H_
#define MATRIXMULT_H_

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/type_traits/is_same.hpp>

#include <geomc/linalg/LinalgTypes.h>
#include <geomc/linalg/mtxdetail/MatrixGlue.h>

namespace geom {
namespace detail {
    
// Mt is a type
#define REQUIRE_MATRIX_T(Mt) \
    typename boost::enable_if< \
        boost::is_base_of< \
            geom::detail::MatrixBase<typename Mt::elem_t, Mt::ROWDIM, Mt::COLDIM, typename Mt::recurring_t>, \
            Mt> \
        >::type

// This is necessary because c++ gets confused by all the different matrix mult 
// implementations which involve permutation matrixes:
#define REQUIRE_NONPERMUTE_MATRIX_T(Mt) \
    typename boost::enable_if_c< \
        boost::is_base_of< \
            geom::detail::MatrixBase<typename Mt::elem_t, Mt::ROWDIM, Mt::COLDIM, typename Mt::recurring_t>, \
            Mt \
        >::value and not \
        boost::is_base_of< \
            geom::PermutationMatrix<Mt::ROWDIM>, \
            Mt \
        >::value \
    >::type


/************************************
 * Matrix mult implementations      *
 ************************************/

// BEWARE: Template ambiguity can pop up and rear its ugly head
//         if you are not careful. Be sure to test these thoroughly.

template <typename Ma, typename Mb, typename Enable=void>
class _ImplMtxMul {
public:
    
    //TODO: this should probably use T = decltype(Ma::elem_t * Mb::elem_t)
    typedef SimpleMatrix<typename Ma::elem_t, Ma::ROWDIM, Mb::COLDIM> return_t;
    
    template <typename Md>
    static void mul(Md *d, const Ma& a, const Mb& b){
        //TODO: There is a trick to optimize *d = *a that convinces the compiler that d != a, thus
        //preventing an expensive fetch operation behind the scenes. We should use this here.
        //google key term would be "aliasing". c++ keyword = "restrict"
        typedef typename Md::elem_t T;
        for (index_t r = 0; r < a.rows(); r++) {
            for (index_t c = 0; c < b.cols(); c++) {
                T sum = 0;
                for (index_t n = 0; n < a.cols(); n++) {
                    sum += a.get(r, n) * b.get(n, c);
                }
                d->set(r, c, sum);
            }
        }
    }
};

template <typename Mat, typename T, index_t M, index_t N>
class _ImplMtxMul<Mat, geom::DiagMatrix<T,M,N>, REQUIRE_MATRIX_T(Mat)> {
public:

    typedef SimpleMatrix<T, Mat::ROWDIM, N> return_t;
    
    // (a x M) * (M x N) -> (a x N)
    template <typename Md>
    static void mul(Md *d, const Mat &a, const geom::DiagMatrix<T,M,N> &b) {
        const T *dp = b.diagonal_begin();
        for (index_t c = 0; c < b.cols(); ++c, ++dp) {
            for (index_t r = 0; r < a.rows(); ++r) {
                d->set(r,c, (*dp) * a.get(r,c));
            }
        }
    }
};

template <typename Mat, typename T, index_t M, index_t N>
class _ImplMtxMul<geom::DiagMatrix<T,M,N>, Mat, REQUIRE_MATRIX_T(Mat)> {
public:

    typedef SimpleMatrix<T, M, Mat::COLDIM> return_t;
    
    // (M x N) * (N x a) -> (N x a)
    template <typename Md>
    static void mul(Md *d, const geom::DiagMatrix<T,M,N> &a, const Mat &b) {
        const T *dp = a.diagonal_begin();
        for (index_t r = 0; r < b.rows(); ++r, ++dp) {
            for (index_t c = 0; c < a.cols(); ++c) {
                d->set(r,c, (*dp) * a.get(r,c));
            }
        }
    }
};

template <typename T, index_t L, index_t M1, index_t M2, index_t N>
class _ImplMtxMul<geom::DiagMatrix<T,L,M1>, geom::DiagMatrix<T,M2,N>, void> {
public:
    
    typedef geom::DiagMatrix<T,L,N> return_t;
    
    template <typename Md>
    static void mul(Md *d, const geom::DiagMatrix<T,L,M1> &a, const geom::DiagMatrix<T,M2,N> &b) {
        T *pa = a.diagonal_begin();
        T *pb = b.diagonal_begin();
        d->setZero();
        //xxx I think this is right? do this on paper.
        for (index_t i = 0; i < a.rows(); ++i, ++pa, ++pb) {
            d->set(i, i, (*pa) * (*pb));
        }
    }
};


// mtx * vec
// (row x col) * (col x 1) = (row x 1)
template <typename Mat, typename T, index_t N>
class _ImplMtxMul<Mat, geom::Vec<T,N>, REQUIRE_MATRIX_T(Mat)> {
public:
    typedef typename _ImplVectorLikeMatrix<T,Mat::ROWDIM,ORIENT_VEC_COL>::type return_t; // Vec<T,M>, unless dynamic.
    
    template <index_t M>
    static void mul(geom::Vec<T,M> *d, const Mat &a, geom::Vec<T,N> b) {
        for (index_t r = 0; r < M; r++) {
            T x = 0;
            for (index_t c = 0; c < N; c++) {
                x += a.get(r,c) * b.get(c);
            }
            (*d)[r] = x;
        }
    }
    
    template <typename Md>
    static void mul(Md *d, const Mat &a, geom::Vec<T,N> b, REQUIRE_MATRIX_T(Md) *dummy=0) {
        for (index_t r = 0; r < a.rows(); r++) {
            T x = 0;
            for (index_t c = 0; c < N; c++) {
                x += a.get(r,c) * b.get(c);
            }
            d->set(r,0,x);
        }
    }
};


// vec * mtx
// (1 x row) * (row x col) = (1 x col)
template <typename Mat, typename T, index_t M>
class _ImplMtxMul<geom::Vec<T,M>, Mat, REQUIRE_MATRIX_T(Mat)> {
public:
    typedef typename _ImplVectorLikeMatrix<T,Mat::COLDIM,ORIENT_VEC_ROW>::type return_t; // Vec<T,N>, unless dynamic
    
    template <index_t N>
    static void mul(geom::Vec<T,N> *d, geom::Vec<T,M> b, const Mat &a) {
        for (index_t c = 0; c < N; c++) {
            T x = 0;
            for (index_t r = 0; r < M; r++) {
                x += a.get(r,c) * b.get(r);
            }
            (*d)[c] = x;
        }
    }
    
    template <typename Md>
    static void mul(Md *d, geom::Vec<T,M> b, const Mat &a, REQUIRE_MATRIX_T(Md) *dummy=0) {
        for (index_t c = 0; c < a.cols(); c++) {
            T x = 0;
            for (index_t r = 0; r < M; r++) {
                x += a.get(r,c) * b.get(r);
            }
            d->set(0,c,x);
        }
    }
};


// mtx * permutation
// (row x col) * (col x col) = (row x col)
// xxx this might fail.
template <typename Mat, index_t N>
class _ImplMtxMul<Mat, geom::PermutationMatrix<N>, REQUIRE_NONPERMUTE_MATRIX_T(Mat)> {
public:
    typedef typename Mat::elem_t T;
    static const index_t M = Mat::ROWDIM;
    
    typedef geom::SimpleMatrix<typename Mat::elem_t,M,N> return_t;
    
    template <typename Md>
    static void mul(Md *d, const Mat &a, const geom::PermutationMatrix<N> &b) {
        index_t rows = a.rows();
        const index_t *p = b.getColSources();
        for (index_t col = 0; col < b.cols(); col++) {
            index_t src_col = p[col];
            std::copy(a.col(src_col), a.col(src_col) + rows, d->col(col));
        }
    }
};


// permutation * mtx
// (row x row) * (row x col) = (row x col)
template <typename Mat, index_t M>
class _ImplMtxMul<geom::PermutationMatrix<M>, Mat, REQUIRE_NONPERMUTE_MATRIX_T(Mat)> {
public:
    typedef typename Mat::elem_t T;
    
    typedef geom::SimpleMatrix<typename Mat::elem_t, M, Mat::COLDIM> return_t;
    
    template <typename Md, index_t N>
    static void mul(Md *d, const geom::PermutationMatrix<N> &a, const Mat &b) {
        index_t cols = b.cols();
        const index_t *p = a.getRowSources();
        for (index_t row = 0; row < a.rows(); row++) {
            index_t src_row = p[row];
            std::copy(b.row(src_row), b.row(src_row) + cols, d->row(row));
        }
    }
};

template <index_t M, index_t N>
class _ImplMtxMul<geom::PermutationMatrix<M>, geom::PermutationMatrix<N>, REQUIRE_MATRIX_T(geom::PermutationMatrix<M>)> { // require is not "neccessary", but helps disambiguate template solving
public:
    
    typedef geom::PermutationMatrix<M> Ma;
    typedef geom::PermutationMatrix<N> Mb;
    typedef geom::PermutationMatrix<((M < N) ? M : N)> return_t;
    
    template <typename Md>
    static void mul(Md *d, const Ma &a, const Ma &b) {
        const index_t *d_a = a.getRowDestinations();
        const index_t *d_b = b.getRowDestinations();
        
        d->setZero();
        for (index_t i = 0; i < a.rows(); i++) {
            d->set(i, d_a[d_b[i]], 1);
        }
    }
    
    template <index_t K>
    static void mul(geom::PermutationMatrix<K> *d, const Ma &a, const Mb &b) {
        const index_t *s_a = a.getRowSources();
        const index_t *s_b = b.getRowSources();
        const index_t *d_a = a.getRowDestinations();
        const index_t *d_b = b.getRowDestinations();
        
        index_t *d_rd = d->getSrcData();
        index_t *d_cd = d->getDstData();
        
        // the best way to understand this is to draw a picture with little 
        // arrows and then look at it.
        for (index_t i = 0; i < a.rows(); i++) {
            d_rd[i] = s_b[s_a[i]];
            d_cd[i] = d_a[d_b[i]];
        }
    }
};

// permutation * vec
// (row x row) * (row x 1) = (row x 1)
template <typename T, index_t N1, index_t N2>
class _ImplMtxMul<geom::PermutationMatrix<N1>, geom::Vec<T,N2>, REQUIRE_MATRIX_T(geom::PermutationMatrix<N1>)> {
public:

    typedef typename _ImplVectorLikeMatrix<T,N2,ORIENT_VEC_COL>::type return_t; // Vec<T,N2>, unless dynamic
    
    static void mul(Vec<T,N2> *d, const geom::PermutationMatrix<N1> &a, geom::Vec<T,N2> b) {
        const index_t *p = a.getRowSources();
        for (index_t i = 0; i < N2; i++) {
            (*d)[i] = b[p[i]];
        }
    }
    
    template <typename Md>
    static void mul(Md *d, const geom::PermutationMatrix<N1> &a, geom::Vec<T,N2> b, REQUIRE_MATRIX_T(Md) *dummy=0) {
        const index_t *p = a.getRowSources();
        for (index_t i = 0; i < N2; i++) {
            d->set(i, 0, b[p[i]]);
        }
    }
};

// vec * permutation
// (1 x col) * (col x col) = (1 x col)
template <typename T, index_t N1, index_t N2>
class _ImplMtxMul<geom::Vec<T,N1>, geom::PermutationMatrix<N2>, REQUIRE_MATRIX_T(geom::PermutationMatrix<N2>)> { // require is not neccessary; this just makes the specialization "more specific"
public:

    typedef typename _ImplVectorLikeMatrix<T,N1,ORIENT_VEC_ROW>::type return_t; // Vec<T,N1>, unless dynamic
    
    static void mul(geom::Vec<T,N1> *d, const geom::Vec<T,N1> &a, const geom::PermutationMatrix<N2> &b) {
        const index_t *p = b.getColSources();
        for (index_t i = 0; i < N1; i++) {
            (*d)[i] = a[p[i]];
        }
    }
    
    template <typename Md>
    static void mul(Md *d, const geom::Vec<T,N1> &a, const geom::PermutationMatrix<N2> &b, REQUIRE_MATRIX_T(Md) *dummy=0) {
        const index_t *p = b.getColSources();
        for (index_t i = 0; i < N1; i++) {
            d->set(0, i, a[p[i]]);
        }
    }
};

//TODO: mult of vector could be made speedy for contiguous matrices.
//      can increment matrix ptr and use     v++ % N.

#undef REQUIRE_MATRIX_T

}; // namespace geom
}; // namespace detail


#endif /* MATRIXMULT_H_ */
