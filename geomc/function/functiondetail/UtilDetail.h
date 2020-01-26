#ifndef _UTILDETAIL_H_
#define _UTILDETAIL_H_

#include <cmath>
#include <type_traits>

namespace geom {
namespace detail {

template <typename T, typename Enable=void>
struct _ImplFMADot {
    static inline T diff_of_products(T a, T b, T c, T d) {
        return a * b - c * d;
    }
    
    static inline T sum_of_products(T a, T b, T c, T d) {
        return a * b + c * d;
    }
};

template <typename T>
struct _ImplFMADot <T, typename std::enable_if<std::is_floating_point<T>::value>::type> {
    static inline T diff_of_products(T a, T b, T c, T d) {
        T cd  = c * d;
        T err = std::fma(-c, d, cd);
        T x   = std::fma(a, b, -cd);
        return x + err;
    }
    
    static inline T sum_of_products(T a, T b, T c, T d) {
        T cd  = c * d;
        T err = std::fma(c, d, -cd);
        T x   = std::fma(a, b,  cd);
        return x + err;
    }
};

}
}

#endif
