#if __cplusplus > 201103L
#include <type_traits>
#endif


namespace geom {


template <typename T, bool Const>
struct ConstType {
    // inaccessible
};

template <typename T>
struct ConstType <T,false> {
    typedef T  type;
    typedef T* pointer_t;
    typedef T& reference_t;
};

template <typename T>
struct ConstType <T,true> {
    typedef const T  type;
    typedef const T* pointer_t;
    typedef const T& reference_t;
};


} // end namespace geom


#if __cplusplus < 201103L

namespace std {
    
    template <bool Cond, typename T, typename F>
    struct conditional {};
    
    template <typename T, typename F>
    struct conditional <true, T, F> {
        typedef T type;
    };
    
    template <typename T, typename F>
    struct conditional <false, T, F> {
        typedef F type;
    };
    
}

#endif