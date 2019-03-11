#include <geomc/linalg/LinalgTypes.h>

// todo: is the StrideIterator faster or slower than the Matrix{Row,Col}Iterators?

namespace geom {

template <typename T, bool Const>
class StrideIterator : 
        public boost::iterator_facade<
            StrideIterator<T,Const>,                               // self type
            T,                                                     // value (pointed to) type
            std::random_access_iterator_tag,                       // implemented concepts (all)
            typename std::conditional<Const,const T&, T&>::type> { // reference to elem type
    
    typedef typename std::conditional<Const,const T*,T*>::type pointer_t;
    
public:
    
    pointer_t base;
    index_t   i;
    index_t   stride;
    index_t   len;
    
    StrideIterator(pointer_t base, index_t i, index_t stride, index_t len):
        base(base),
        i(i),
        stride(stride),
        len(len) {}
    
    template <bool K>
    inline bool equal(const StrideIterator<T,K>& other) const {
        return other.base == base and other.i == i;
    }
    
    inline bool equal(const T* other) const {
        return base + i == other;
    }
    
    inline void increment() {
        // we don't call advance(1) because handling the negative offset case
        // requires a few extra operations (including a branch), and this is
        // very inner-loop code.
        index_t i_1 = i + stride;
        i = i_1 % len + i_1 / len;
    }
    
    inline void decrement() {
        advance(-1);
    }
    
    inline void advance(index_t dx) {
        index_t i_1 = i + stride * dx;
        i = i_1 % len + i_1 / len; // needs floor(i/len), not truncate(i/len).
        if (i_1 < 0) {
            // C and C++ integer division truncate towards zero 
            // instead of rounding to negative infinity, hence the `-1`.
            // (this is arguably a deep design flaw in C and C++. 
            // Guido explains why Python took the high road:
            // http://python-history.blogspot.com/2010/08/why-pythons-integer-division-floors.html)
            // Is there a branchless / caseless way to handle this? I'd love to know.
            i += len - 1;
        }
    }
    
    inline typename std::conditional<Const,const T&, T&>::type dereference() const {
        return *(base + i);
    }
    
    template <bool K>
    inline index_t distance_to(const StrideIterator<T,K>& other) const {
        return (other.i % stride - i % stride) * (len / stride) + 
               (other.i / stride - i / stride);
    }
    
};


template <typename T, MatrixLayout L>
class FlatMatrixLayout {
public:
    
    typedef                const T* const_row_iterator;
    typedef  StrideIterator<T,true> const_col_iterator;
    typedef                      T* row_iterator;
    typedef StrideIterator<T,false> col_iterator;
    
    static inline index_t index(index_t r, index_t c, index_t rows, index_t cols) {
        return r * cols + c;
    }
    
    template <bool Const>
    static inline typename std::conditional<Const, const T*, T*>::type row(
            typename std::conditional<Const, const T*, T*>::type base, index_t r, index_t rows, index_t cols) {
        return base + r * cols;
    }
    
    template <bool Const>
    static inline StrideIterator<T,Const> col(
            typename std::conditional<Const, const T*, T*>::type base, index_t c, index_t rows, index_t cols) {
        return StrideIterator<T,Const>(base, c, cols, rows * cols);
    }
    
};


template <typename T>
class FlatMatrixLayout<T, COL_MAJOR> {
public:
    
    typedef  StrideIterator<T,true> const_row_iterator;
    typedef                const T* const_col_iterator;
    typedef StrideIterator<T,false> row_iterator;
    typedef                      T* col_iterator;
    
    static inline index_t index(index_t r, index_t c, index_t rows, index_t cols) {
        return c * rows + r;
    }
    
    template <bool Const>
    static inline typename std::conditional<Const, const T*, T*>::type col(
            typename std::conditional<Const, const T*, T*>::type base, index_t c, index_t rows, index_t cols) {
        return base + c * rows;
    }
    
    template <bool Const>
    static inline StrideIterator<T,Const> row(
            typename std::conditional<Const, const T*, T*>::type base, index_t r, index_t rows, index_t cols) {
        return StrideIterator<T,Const>(base, r, rows, rows * cols);
    }
    
};

};