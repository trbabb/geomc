#include <geomc/linalg/LinalgTypes.h>

namespace geom {


template <typename T, bool Const>
struct TransposeIterator :
        public boost::iterator_facade
        <
            TransposeIterator<T,Const>,      // self_t
            T,                               // pointed-to type
            std::random_access_iterator_tag, // implemented concepts
                                             // ↓ elem type
            typename std::conditional<Const, const T&, T&>::type
        >
{
    
    typedef typename std::conditional<Const, const T*, T*>::type pointer_t;
    typedef typename std::conditional<Const, const T&, T&>::type ref_t;
    
    /*
        maj--->         min
         0  1  2  3  4  |
         5  6  7  8  9  v
        10 11 12 13 14
        [] <-- "end"
    */
    
    const pointer_t base;
    const index_t   major; // # elements along major (consecutive) axis
    const index_t   minor; // # elements along minor axis
          index_t   i;     // ordinal of this iterator
          pointer_t p;     // pointed memory location
    
    TransposeIterator(pointer_t base, index_t major, index_t minor, index_t i=0):
        base(base),
        major(major),
        minor(minor),
        i(i),
        p(offs()) {}
    
    // compute the new memory location
    inline pointer_t offs() const {
        // handle "off-end" correctly:
        index_t s = major * minor;
        index_t j = i % s;
        index_t k = i / s;
        // compute transposed offset:
        return base + (j % minor) * major + (j / minor) + k * s;
    }
    
    inline bool equal(const T* other) const {
        return this->p == other;
    }
    
    template <bool K>
    inline bool equal(const TransposeIterator<T,K>& other) const {
        return this->p == other.p;
    }
    
    inline void increment() {
        this->i += 1;
        this->p  = offs();
    }
    
    inline void decrement() {
        this->i -= 1;
        this->p  = offs();
    }
    
    inline void advance(index_t dx) {
        this->i += dx;
        this->p  = offs();
    }
    
    inline ref_t dereference() const {
        return *p;
    }
    
    template <bool K>
    inline index_t distance_to(const TransposeIterator<T,K>& other) const {
        index_t   j = other.p - other.base;
        index_t o_i = (j % major) * minor + (j / major);
        
        return o_i - i;
    }
    
};


template <typename T, MatrixLayout L>
class FlatMatrixLayout {
public:
    
    typedef                   const T* const_row_iterator;
    typedef TransposeIterator<T, true> const_col_iterator;
    typedef                         T* row_iterator;
    typedef TransposeIterator<T,false> col_iterator;
    
    static inline index_t index(index_t r, index_t c, index_t rows, index_t cols) {
        return r * cols + c;
    }
    
    template <bool Const>
    static inline typename std::conditional<Const, const T*, T*>::type 
        row(
            typename std::conditional<Const, const T*, T*>::type base, 
            index_t r,
            index_t rows,
            index_t cols)
    {
        return base + r * cols;
    }
    
    template <bool Const>
    static inline TransposeIterator<T,Const> col(
            typename std::conditional<Const, const T*, T*>::type base, 
            index_t c, 
            index_t rows, 
            index_t cols)
    
    {
        return TransposeIterator<T,Const>(base, cols, rows, c * rows);
    }
    
};


template <typename T>
class FlatMatrixLayout<T, COL_MAJOR> {
public:
    
    typedef TransposeIterator<T, true> const_row_iterator;
    typedef                   const T* const_col_iterator;
    typedef TransposeIterator<T,false> row_iterator;
    typedef                         T* col_iterator;
    
    static inline index_t index(index_t r, index_t c, index_t rows, index_t cols) {
        return c * rows + r;
    }
    
    template <bool Const>
    static inline typename std::conditional<Const, const T*, T*>::type 
        col(
            typename std::conditional<Const, const T*, T*>::type base,
            index_t c,
            index_t rows,
            index_t cols)
    {
        return base + c * rows;
    }
    
    template <bool Const>
    static inline TransposeIterator<T,Const> 
        row(
            typename std::conditional<Const, const T*, T*>::type base,
            index_t r,
            index_t rows,
            index_t cols)
    {
        return TransposeIterator<T,Const>(base, rows, cols, r * cols);
    }
    
};

} // namespace geom
