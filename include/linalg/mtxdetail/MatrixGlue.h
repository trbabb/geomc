/*
 * MatrixGlue.h
 *
 *  Created on: Jan 28, 2013
 *      Author: tbabb
 */

#ifndef MATRIXGLUE_H_
#define MATRIXGLUE_H_

#include <boost/smart_ptr/shared_array.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include "shape/GridIterator.h"
#include "MatrixMacro.h"

namespace geom {

typedef Rect<index_t,2> MatrixRegion;
typedef Vec<index_t,2>  MatrixCoord;

namespace detail {

/****************************************************
 * Storage class                                    *
 *                                                  *
 * Matrices may be of either static or dynamic      *
 * dimensions. Statically-sized matrices may use    *
 * statically-allocated storage and avoid expensive *
 * mallocs() upon creation and duplication. This    *
 * class allows matrix templates to choose between  *
 * the two methods of storage based on <N>.         *
 ****************************************************/

template <typename T, index_t N>
struct Storage {
    T data[N];
    
    explicit Storage(index_t n) {}
    
    inline       T* get()       { return data; }
    inline const T* get() const { return data; }
    inline index_t size() const { return N; }
};

template <typename T>
struct Storage<T,DYNAMIC_DIM> {
    boost::shared_array<T> data;
    
    explicit Storage(index_t n):
                    data(new T[n]) {}
    
    inline       T* get()       { return data.get(); }
    inline const T* get() const { return data.get(); }
};

/****************************************************
 * Temporary storage class                          *
 *                                                  *
 * Like Storage, except uses bare array pointers    *
 * internally with no reference counting. Data will *
 * be deallocated automatically on destruction. For *
 * use when ownership never changes hands.          *
 ****************************************************/

template <typename T, index_t N>
struct TemporaryStorage {
    T data[N];
    
    explicit TemporaryStorage(index_t n) {}
    
    inline       T* get()       { return data; }
    inline const T* get() const { return data; }
    inline index_t size() const { return N; }
    
    inline T& operator[](index_t idx)       { return data[idx]; }
    inline T  operator[](index_t idx) const { return data[idx]; }
};

template <typename T>
struct TemporaryStorage<T,DYNAMIC_DIM> {
    T *data;
    index_t sz;
    
    explicit TemporaryStorage(index_t n):
        data(new T[n]),
        sz(n) {}
    
    ~TemporaryStorage() {
        delete [] data;
    }
    
    inline       T* get()       { return data; }
    inline const T* get() const { return data; }
    inline index_t size() const { return sz; }
    
    inline T& operator[](index_t idx)       { return data[idx]; }
    inline T  operator[](index_t idx) const { return data[idx]; }
}; 

/****************************************************
 * Dimension classes                                *
 *                                                  *
 * Let's say we have a Storage object, which may    *
 * wrap either a static or a dynamic array. If the  *
 * array is dynamic, we must store the length. We   *
 * don't want to pay for that length counter if the *
 * wrapped memory is a static array, however, since *
 * the length is known at compile time. This class  *
 * allocates storage for one size counter if <N> is *
 * DYNAMIC_DIMENSION, otherwise no storage is       *
 * allocated and <N> is the reported array length.  *
 ****************************************************/

template <index_t N>
struct Dimension {
    typedef index_t storage_t[0];
    static inline void set(storage_t s, index_t v) {}
    static inline index_t value(const storage_t v) { return N; }
    
};

template <>
struct Dimension<DYNAMIC_DIM> {
    typedef index_t storage_t[1];
    static inline void set(storage_t s, index_t v) { s[0] = v; }
    static inline index_t value(const storage_t v) { return *v; }
    
};

/****************************************************
 * Conditional template helper                      *
 *                                                  *
 * DefinedIf<cond,val>::value is only defined if    *
 * cond is true, and takes on the value <val>.      *
 *                                                  *
 * Used by SimpleMatrix and similar classes to init *
 * their default constructor args with dimensions.  *   
 * (If a dimension is dynamic, we should not allow  *
 * a default length for that dimension, and this    *
 * template disables said default).                 *
 ****************************************************/

template <bool cond, index_t val=0>
struct DefinedIf {
    typedef void type;
};

template <index_t val>
struct DefinedIf<true, val> {
    typedef index_t type;
    
    static const index_t value = val;
};

/******************************************************
 * Storage object count templates                     *
 *                                                    *
 * These are used to track whether a matrix has a     *
 * static or dynamic storage object count, which in   *
 * turn affects the type of the storage ID buffer     *
 * (used for tracking memory aliasing).               *
 *                                                    *
 * We have a special type for the storage buffer      *
 * because we would prefer not to malloc() any time   *
 * we check for memory aliasing (i.e. every matrix    *
 * mult) unless we absolutely have to. Static logic   *
 * about buffer size allows for static array allocs.  *
 *                                                    *
 * A compound matrix, like AugmentedMatrix, for       *
 * example, may have either static or dynamic storage *
 * count.                                             *
 ******************************************************/

template <typename M>
struct _ImplStorageObjCount {
    static const index_t count = 1;  // to be overridden by mtx-specific specializations
};

/******************************************************
 * IsMatrix check                                     *
 ******************************************************/

template <typename Mt, typename Enable=void>
struct IsMatrix {
    static const bool val = false;
};

template <typename Mt>
struct IsMatrix<Mt, typename boost::enable_if_c<
            boost::is_base_of<
                geom::detail::MatrixBase<typename Mt::elem_t, Mt::ROWDIM, Mt::COLDIM, Mt>, 
                Mt>::value,
            void>::type> {
    static const bool val = true;
};

/*****************************************************
 * Matrix dimension agreement                        *
 *                                                   *
 * Check if two matrices' dimensions match           *
 *****************************************************/

template <typename Ma, typename Mb, typename Enable=void>
struct MatrixDimensionMatch {
    static const bool isStaticMatch = false;
    
    static inline bool isMatch(const Ma &ma, const Mb &mb) {
        return false;
    }
    
    static inline void check(const Ma &ma, const Mb &mb) {
        throw DimensionMismatchException(ma.rows(), ma.cols(), mb.rows(), mb.cols());
    }
};

template <typename Ma, typename Mb>
struct MatrixDimensionMatch <Ma, Mb, 
        typename boost::enable_if_c<
            ((Ma::ROWDIM == Mb::ROWDIM or Ma::ROWDIM * Mb::ROWDIM == 0) and 
             (Ma::COLDIM == Mb::COLDIM or Ma::COLDIM * Mb::COLDIM == 0))
        >::type> {
    
    static const bool isStaticMatch = true;
    
    static inline bool isMatch(const Ma &ma, const Mb &mb) {
        if ((Ma::ROWDIM * Mb::ROWDIM == 0 and ma.rows() != mb.rows()) or
            (Ma::COLDIM * Mb::COLDIM == 0 and ma.cols() != mb.cols())) {
            return false;
        }
        return true;
    }
    
    static inline void check(const Ma& ma, const Mb &mb) {
        if (not isMatch(ma, mb)) {
            throw DimensionMismatchException(ma.rows(), ma.cols(), mb.rows(), mb.cols());
        }
    }
};

/******************************************************
 * Vector-like object instantiator                    *
 *                                                    *
 * For all matrix-related functions that can produce  *
 * a Vec, we must also allow for the case where all   *
 * dimensions are dynamic. In this case, instead of a *
 * Vec we return a 1 x N dynamic SimpleMatrix.        *
 * If we ever implement dynamic-length vectors, this  *
 * will become unnecessary.                           *
 ******************************************************/

template <typename T, index_t N, VecOrientation O>
struct _ImplVectorLikeMatrix {
    typedef geom::Vec<T,N> type;
};

template <typename T>
struct _ImplVectorLikeMatrix<T, DYNAMIC_DIM, ORIENT_VEC_ROW> {
    typedef geom::SimpleMatrix<T, 1, DYNAMIC_DIM> type;
};

template <typename T>
struct _ImplVectorLikeMatrix<T, DYNAMIC_DIM, ORIENT_VEC_COL> {
    typedef geom::SimpleMatrix<T, DYNAMIC_DIM, 1> type;
};

template <typename T>
struct _ImplVectorLikeMatrix<T, DYNAMIC_DIM, ORIENT_VEC_UNKNOWN> {
    typedef geom::SimpleMatrix<T, DYNAMIC_DIM, DYNAMIC_DIM> type;
};

/******************************************************
 * Dimension checking adaptors                        *
 *                                                    *
 * Templates to query the row/column dimensions of a  *
 * matrix-or-vector. Allows matrices/vectors to be    *
 * interoperable. Since vectors do not inherently     *
 * know whether they are row or column matrices, the  *
 * choice is made by through the adaptor.             *
 *                                                    *
 * These also have the function of verifying that     *
 * some object is a linear algebra object.            *
 *                                                    *
 * Primarily used by mul(mtx,mtx).                    *
 ******************************************************/

template <typename T, VecOrientation O, typename Enable=void>
struct _ImplMtxAdaptor {
    // pass
};

template <typename Mx, VecOrientation O>
struct _ImplMtxAdaptor <Mx, O, typename boost::enable_if_c<IsMatrix<Mx>::val, void>::type> {
    static const index_t ROWDIM = Mx::ROWDIM;
    static const index_t COLDIM = Mx::COLDIM;
    
    static index_t rows(const Mx &m) { return m.rows(); }
    static index_t cols(const Mx &m) { return m.cols(); }
};

template <typename T, index_t N>
struct _ImplMtxAdaptor <geom::Vec<T,N>, ORIENT_VEC_COL, void> {
    static const index_t ROWDIM = N;
    static const index_t COLDIM = 1;
    
    static index_t rows(const Vec<T,N>& v) { return ROWDIM; }
    static index_t cols(const Vec<T,N>& v) { return COLDIM; }
};

template <typename T, index_t N>
struct _ImplMtxAdaptor <geom::Vec<T,N>, ORIENT_VEC_ROW, void> {
    static const index_t ROWDIM = 1;
    static const index_t COLDIM = N;
    
    static index_t rows(const Vec<T,N>& v) { return ROWDIM; }
    static index_t cols(const Vec<T,N>& v) { return COLDIM; }
};

template <typename T, index_t N>
struct _ImplMtxAdaptor <geom::Vec<T,N>, ORIENT_VEC_UNKNOWN, void> {
    static const index_t ROWDIM = DYNAMIC_DIM;
    static const index_t COLDIM = DYNAMIC_DIM;
};

///////// orientation finder /////////
// the result of a matrix * vector or vector * matrix
// operation may be either a row or a column vector.
// this adaptor deciphers the orientation based
// on the static dimensionality of the operands.

template <typename Ma, typename Mb, typename Md, typename Enable=void>
struct _ImplVecMulOrient {
    static const VecOrientation orient = ORIENT_VEC_UNKNOWN;
};

template <typename Ma, typename Mb, typename Md>
struct _ImplVecMulOrient <Ma, Mb, Md, typename boost::enable_if_c<
            _ImplMtxAdaptor<Ma, ORIENT_VEC_ROW>::ROWDIM == _ImplMtxAdaptor<Md, ORIENT_VEC_COL>::ROWDIM and 
            _ImplMtxAdaptor<Mb, ORIENT_VEC_COL>::COLDIM == _ImplMtxAdaptor<Md, ORIENT_VEC_COL>::COLDIM, 
            void>::type> {
    static const VecOrientation orient = ORIENT_VEC_COL;
};

template <typename Ma, typename Mb, typename Md>
struct _ImplVecMulOrient <Ma, Mb, Md, typename boost::enable_if_c<
            _ImplMtxAdaptor<Ma, ORIENT_VEC_ROW>::ROWDIM == _ImplMtxAdaptor<Md, ORIENT_VEC_ROW>::ROWDIM and 
            _ImplMtxAdaptor<Mb, ORIENT_VEC_COL>::COLDIM == _ImplMtxAdaptor<Md, ORIENT_VEC_ROW>::COLDIM and not
            /* avoid ambiguity with column orienation template: */
           (_ImplMtxAdaptor<Ma, ORIENT_VEC_ROW>::ROWDIM == _ImplMtxAdaptor<Md, ORIENT_VEC_COL>::ROWDIM and 
            _ImplMtxAdaptor<Mb, ORIENT_VEC_COL>::COLDIM == _ImplMtxAdaptor<Md, ORIENT_VEC_COL>::COLDIM), 
            void>::type> {
    static const VecOrientation orient = ORIENT_VEC_ROW;
};

///////// orientation finder /////////
// how shall we orient two possible-vectors such that
// their static dimensions agree?
// if the orientation cannot be proven (and is not irrelevant,
// as in the case of two vectors, or two matrices), then
// ORIENT_VEC_UNKNOWN is resolved.

template <typename Ma, typename Mb, typename Enable=void>
struct _ImplVecOrient {
    static const VecOrientation orient = ORIENT_VEC_UNKNOWN;
};

// vec(N x 1) == mtx(N x 1)
// col vector and matrix
template <typename T, index_t N, typename Mb>
struct _ImplVecOrient <Vec<T,N>, Mb, typename boost::enable_if_c<
                (N == Mb::ROWDIM or Mb::ROWDIM == DYNAMIC_DIM) and
                (Mb::COLDIM == 1 or Mb::COLDIM == DYNAMIC_DIM) and not
                (Mb::COLDIM == DYNAMIC_DIM and Mb::ROWDIM == DYNAMIC_DIM),
            void>::type> {
    static const VecOrientation orient = ORIENT_VEC_COL;
};

// mtx(N x 1) == vec(N x 1)
// matrix and col vector
template <typename T, index_t N, typename Mb>
struct _ImplVecOrient <Mb, Vec<T,N>, typename boost::enable_if_c<
                (N == Mb::ROWDIM or Mb::ROWDIM == DYNAMIC_DIM) and
                (Mb::COLDIM == 1 or Mb::COLDIM == DYNAMIC_DIM) and not
                (Mb::COLDIM == DYNAMIC_DIM and Mb::ROWDIM == DYNAMIC_DIM),
            void>::type> {
    static const VecOrientation orient = ORIENT_VEC_COL;
};

// vec(1 x N) == mtx(1 x N)
// row vector and matrix
template <typename T, index_t N, typename Mb>
struct _ImplVecOrient <Vec<T,N>, Mb, typename boost::enable_if_c<
                (N == Mb::COLDIM or Mb::COLDIM == DYNAMIC_DIM) and
                (Mb::ROWDIM == 1 or Mb::ROWDIM == DYNAMIC_DIM) and not
                (Mb::COLDIM == DYNAMIC_DIM and Mb::ROWDIM == DYNAMIC_DIM),
            void>::type> {
    static const VecOrientation orient = ORIENT_VEC_ROW;
};

// mtx(1 x N) == vec(1 x N)
// matrix and row vector
template <typename T, index_t N, typename Mb>
struct _ImplVecOrient <Mb, Vec<T,N>, typename boost::enable_if_c<
                (N == Mb::COLDIM or Mb::COLDIM == DYNAMIC_DIM) and
                (Mb::ROWDIM == 1 or Mb::ROWDIM == DYNAMIC_DIM) and not
                (Mb::COLDIM == DYNAMIC_DIM and Mb::ROWDIM == DYNAMIC_DIM),
            void>::type> {
    static const VecOrientation orient = ORIENT_VEC_ROW;
};

// two vectors
template <typename A, typename B, index_t N>
struct _ImplVecOrient <Vec<A,N>, Vec<B,N>, void> {
    // arbitrary
    static const VecOrientation orient = ORIENT_VEC_COL;
};

// two matrices
template <typename Ma, typename Mb>
struct _ImplVecOrient <Ma, Mb, typename boost::enable_if_c<(IsMatrix<Ma>::val and IsMatrix<Mb>::val), void>::type> {
    // arbitrary
    static const VecOrientation orient = ORIENT_VEC_COL;
};

/////////////////////////////////////////
// Are two matrices multipliable?      //
/////////////////////////////////////////

template <typename Ma, typename Mb, typename Enable=void> 
struct MatrixMultipliable {
    static const bool val = false;
};

template <typename Ma, typename Mb>
struct MatrixMultipliable <Ma, Mb, 
    typename boost::enable_if_c<
            /* inner dimension match */
            (_ImplMtxAdaptor<Ma, ORIENT_VEC_ROW>::COLDIM == _ImplMtxAdaptor<Mb, ORIENT_VEC_COL>::ROWDIM or
            /* ...or dynamic inner dimension demands runtime check: */
            _ImplMtxAdaptor<Ma, ORIENT_VEC_ROW>::COLDIM == DYNAMIC_DIM or
            _ImplMtxAdaptor<Mb, ORIENT_VEC_COL>::ROWDIM == DYNAMIC_DIM),
//MATRIX_MUL_DIM_AGREE(Ma,Mb),
        void>::type> {
    static const bool val = true;
};

/////////////////////////////////////////
// Is a matrix a valid destination for //
// a matrix mult operation?            //
/////////////////////////////////////////

#define MATRIX_DEST_MATCH(Ma,Mb,Md) \
    /* result row dimension matches: */ \
   ((detail::_ImplMtxAdaptor<Md, detail::_ImplVecMulOrient<Ma,Mb,Md>::orient>::ROWDIM == detail::_ImplMtxAdaptor<Ma, detail::ORIENT_VEC_ROW>::ROWDIM or \
    /* ...or result row dimension is dynamic: */ \
     detail::_ImplMtxAdaptor<Md, detail::_ImplVecMulOrient<Ma,Mb,Md>::orient>::ROWDIM  * detail::_ImplMtxAdaptor<Ma, detail::ORIENT_VEC_ROW>::ROWDIM == 0) and \
    /* column dimension matches: */ \
    (detail::_ImplMtxAdaptor<Md, detail::_ImplVecMulOrient<Ma,Mb,Md>::orient>::COLDIM == detail::_ImplMtxAdaptor<Mb, detail::ORIENT_VEC_COL>::COLDIM or \
    /* or result column dimension is dynamic: */ \
     detail::_ImplMtxAdaptor<Md, detail::_ImplVecMulOrient<Ma,Mb,Md>::orient>::COLDIM  * detail::_ImplMtxAdaptor<Mb, detail::ORIENT_VEC_COL>::COLDIM == 0))

template <typename Ma, typename Mb, typename Md, typename Enable=void>
struct _ImplMtxResult {
    static const MatrixResultAgreement agreement = MTX_RESULT_MISMATCH;
};

template <typename Ma, typename Mb, typename Md>
struct _ImplMtxResult <Ma, Mb, Md,
       typename boost::enable_if_c<MATRIX_DEST_MATCH(Ma,Mb,Md) and detail::_ImplVecMulOrient<Ma,Mb,Md>::orient != ORIENT_VEC_UNKNOWN, void>::type> {
    static const MatrixResultAgreement agreement = MTX_RESULT_MATCH;
};

template <typename Ma, typename Mb, typename Md>
struct _ImplMtxResult <Ma, Mb, Md,
       typename boost::enable_if_c<MATRIX_DEST_MATCH(Ma,Mb,Md) and detail::_ImplVecMulOrient<Ma,Mb,Md>::orient == ORIENT_VEC_UNKNOWN, void>::type> {
   static const MatrixResultAgreement agreement = MTX_RESULT_UNKNOWN;
};

#undef MATRIX_DEST_MATCH

}; /* namespace detail */
}; /* namespace geom */



#endif /* MATRIXGLUE_H_ */
