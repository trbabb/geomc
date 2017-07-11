/*
 * MatrixIterator.h
 *
 *  Created on: Jun 25, 2013
 *      Author: tbabb
 */

#ifndef MATRIXITERATOR_H_
#define MATRIXITERATOR_H_

namespace geom {
namespace detail {

/****************************************************
 * SubsetIterator class                             *
 *                                                  *
 * Used for iterating over sub-regions of matrices. *
 * Two SubsetIterators are equal if they refer to   *
 * the same element of the same matrix. Iteration   *
 * happens in row-major order.                      *
 ****************************************************/

//PONDER: can we factor common stuff out somehow?
//        > consider a unified matrix template which implements
//          its functionality by operating with iterators.
//          it could be templated over the type of the iterator.
//          consider how this may allow row/column major memory layout
//          (mtx returns separate iterators row/col iteration, some of which
//          may be native ptrs, e.g.)

//TODO: diagonal iterators.
//TODO: nonzero iterators?
//TODO: verify behavior of i->x()
//TODO: make a special randaccess iterator for contiguous types, such that i[x] is a bit faster? 
//      i.e. use a ptr rather than a point.


// PROBLEM: We are looping in the wrong order.
//    I believe because row, col is a funny order.

template <typename M, typename RefType>
class MtxSubsetIterator : public boost::iterator_facade<MtxSubsetIterator<M,RefType>, // self type
                                             typename M::elem_t,                      // value (pointed to) type
                                             std::random_access_iterator_tag,         // implemented concepts (all)
                                             RefType> {                               // reference to elem type (may be a proxy)
    
public:
    
    // TODO: optimize this class. maybe counter math is faster?
    
    typedef RefType ref_t;

    M *mtx;
    GridIterator<index_t,2,ARRAYORDER_LAST_DIM_CONSECUTIVE> p;
    
    MtxSubsetIterator(M *m, const MatrixRegion &region):
                          mtx(m),
                          p(region) {}
    
    MtxSubsetIterator(M *m, const MatrixRegion &region, const MatrixCoord &pt):
                          mtx(m),
                          p(region, pt) {}
    
    MtxSubsetIterator(M *m, MatrixCoord pt=MatrixCoord::zeros):
                          mtx(m),
                          p(MatrixRegion(MatrixCoord::zeros,
                                         MatrixCoord(m->rows(), m->cols())),
                            pt) {}
    
    MtxSubsetIterator(M *m, const GridIterator<index_t,2,ARRAYORDER_LAST_DIM_CONSECUTIVE> &p):
                          mtx(m),
                          p(p) {}
public:
    
    MatrixCoord point() {
        return *p;
    }
    
    inline MtxSubsetIterator<M,RefType> begin() {
        return MtxSubsetIterator(mtx, p.begin());
    }
    
    inline MtxSubsetIterator<M,RefType> end() {
        return MtxSubsetIterator(mtx, p.end()); 
    }
    
    
private:
    
    friend class boost::iterator_core_access;
    
    inline bool equal(const MtxSubsetIterator<M,RefType> &other) const {
        return (other.mtx == mtx) && (other.p == p);
    }
    
    inline void increment() { 
        p++;
    }
    
    inline void decrement() {
        p--;
    }
    
    inline void advance(index_t n) {
        p += n;
    }
    
    inline index_t distance_to(const MtxSubsetIterator<M,RefType> &other) const {
        return other.p - p; //this right?
    }
    
    inline RefType dereference() const {
        MatrixCoord pt = *p;
        return mtx->get(pt.row, pt.col); //make sure this is calling the function with the right const-ness!
    }
    
};

/*****************************************************
 * Col iterator                                      *
 *****************************************************/

template <typename M, typename RefType>
class MtxColIterator : public boost::iterator_facade<MtxColIterator<M,RefType>, // self type
                                             typename M::elem_t,                // value (pointed to) type
                                             std::random_access_iterator_tag,   // implemented concepts (all)
                                             RefType> {                         // reference to elem type (may be a proxy)
    
public:
    
    typedef RefType ref_t;

    M *mtx;
    MatrixCoord pt;
    
    MtxColIterator(M *m, const MatrixCoord &pt):
                         mtx(m),
                         pt(pt) {}
    MtxColIterator(M *m, index_t col):
                         mtx(m),
                         pt(0,col) {}
    
private:
    
    friend class boost::iterator_core_access;
    
    inline bool equal(const MtxColIterator<M,RefType> &other) const {
        return (other.mtx == mtx) && (other.pt == pt);
    }
    
    inline void increment() {
        pt.row++;
    }
    
    inline void decrement() {
        pt.row--;
    }
    
    inline void advance(index_t n) {
        pt.row += n;
    }
    
    inline index_t distance_to(const MtxColIterator<M,RefType> &other) const {
        return other.pt.row - pt.row;
    }
    
    inline RefType dereference() const {
        return mtx->get(pt.row, pt.col); //make sure this is calling the function with the right const-ness!
    }
    
};

/*****************************************************
 * Row iterator                                      *
 *****************************************************/

template <typename M, typename RefType>
class MtxRowIterator : public boost::iterator_facade<MtxRowIterator<M,RefType>,       // self type
                                                  typename M::elem_t,              // value (pointed to) type
                                                  std::random_access_iterator_tag, // implemented concepts (all)
                                                  RefType> {                       // reference to elem type (may be a proxy)
    
public:
    
    typedef RefType ref_t;

    M *mtx;
    MatrixCoord pt;
    
    MtxRowIterator(M *m, const MatrixCoord &pt):
                         mtx(m),
                         pt(pt) {}
    MtxRowIterator(M *m, index_t row):
                         mtx(m),
                         pt(row,0) {}
    
private:
    
    friend class boost::iterator_core_access;
    
    inline bool equal(const MtxColIterator<M,RefType> &other) const {
        return (other.mtx == mtx) && (other.pt == pt);
    }
    
    inline void increment() {
        pt.col++;
    }
    
    inline void decrement() {
        pt.col--;
    }
    
    inline void advance(index_t n) {
        pt.col += n;
    }
    
    inline index_t distance_to(const MtxRowIterator<M,RefType> &other) const {
        return other.pt.col - pt.col;
    }
    
    inline RefType dereference() const {
        return mtx->get(pt.row, pt.col); // always make sure this is calling the function with the right const-ness!
    }
    
};

/*****************************************************
 * AssignmentProxy class                             *
 *                                                   *
 * Returned by certain matrix methods in lieu of a   *
 * memory reference, when a memory reference doesn't *
 * make sense. For example, by double-indexing       *
 * a sparse matrix (a memory location for the        *
 * requested element may not exist yet), or by       *
 * dereferencing iterators provided by such classes. *
 * Assigning to these objects inserts an element     *
 * into the associated matrix, if possible.          *
 *****************************************************/

template <typename M, typename T> class MtxAssignmentProxy {
    M *m;
    index_t row, col;
public:
    
    MtxAssignmentProxy(M *m, index_t row, index_t col):m(m), row(row), col(col) {}
    
    // copy assignment
    // like a reference, assign through to the pointed object.
    // contrary to convention, we do NOT mirror this action in the copy constructor
    // this is to mirror the behavior of references, whose ptr can only be assigned at construction.
    inline MtxAssignmentProxy<M,T>& operator=(const MtxAssignmentProxy<M,T> &other) {
        if (&other != this) {
            m->set(row, col, other);
        }
        return *this;
    }
    
    inline const typename M::elem_t operator=(typename M::elem_t val) { //return const T& for some mtxs, T for those which use temporaries?
        m->set(row, col, val);
        return val;
    }
    
    //TODO: this should be a reference?
    inline operator T() const {
        return ((const M*)m)->get(row, col);
    }
    
};

/******************************************************
 * Iterator type traits classes                       *
 *                                                    *
 * Used to inform the matrix base class about the     *
 * derived class's choice of iterator implementation. *
 *                                                    *
 * Because c++ is a stupid, broken language, it       *
 * cannot resolve type declaration dependencies if    *
 * the type definitions belong to the classes         *
 * themselves.                                        *
 *   http://tinyurl.com/cvtlqow                       *
 *   http://tinyurl.com/bl8aznk                       *
 *                                                    *
 * C++ cannot reason about types that belong to other *
 * "incomplete" types, even if a complete definition  *
 * for that sub-type already exists and is itself     *
 * fully defined. If instead C++ built a graph of     *
 * type declarations and topologically sorted their   *
 * dependencies, the problem would be trivially       *
 * solved. But alas, attacking problems in such a     *
 * sensible way is not part of the c++ design         * 
 * philosophy.                                        *
 ******************************************************/

// we have to pass T (rather than use M::elem_t) because, again, 
// c++ is a special child and can't type-solve its way out of wet paper bag.
// specifically, when these are instantiated, M is an incomplete type.
template <typename M, typename T>
struct _ImplMtxReftype {
    typedef MtxAssignmentProxy<M, T> reference;
};

template <typename M, typename T>
struct _ImplMtxConstRowtype {
    typedef MtxRowIterator<M, const T> const_row_iterator;
};

template <typename M, typename T>
struct _ImplMtxRowtype {
    private:
    typedef typename _ImplMtxReftype<M, T>::reference impl_reference;
    
    public:
    typedef MtxRowIterator<M,impl_reference> row_iterator;
};

}; // namespace detail
}; // namespace geom


#endif /* MATRIXITERATOR_H_ */
