#ifndef _UTILDETAIL_H_
#define _UTILDETAIL_H_

#include <cmath>
#include <type_traits>
#include <geomc/function/Dual.h>

namespace geom {
namespace detail {

template <typename T, typename Enable=void>
struct _ImplFMADot {
    static constexpr T diff_of_products(T a, T b, T c, T d) {
        return a * b - c * d;
    }
    
    static constexpr T sum_of_products(T a, T b, T c, T d) {
        return a * b + c * d;
    }
    
    static constexpr T multiply_add(T a, T b, T c) {
        return a * b + c;
    }
};

template <typename T>
struct _ImplFMADot <T, typename std::enable_if<std::is_floating_point<T>::value>::type> {
    static inline T diff_of_products(T a, T b, T c, T d) {
        T cd  = c * d;
        T err = std::fma(-c, d,  cd);
        T x   = std::fma( a, b, -cd);
        return x + err;
    }
    
    static inline T sum_of_products(T a, T b, T c, T d) {
        T cd  = c * d;
        T err = std::fma(c, d, -cd);
        T x   = std::fma(a, b,  cd);
        return x + err;
    }
    
    static inline T multiply_add(T a, T b, T c) {
        return std::fma(a, b, c);
    }
};

template <typename T, DiscontinuityPolicy P>
struct _ImplFMADot<Dual<T,P>, void> {
    static inline Dual<T,P> diff_of_products(Dual<T,P> a, Dual<T,P> b, Dual<T,P> c, Dual<T,P> d) {
        T dp = sum_of_products(a.x, b.dx, a.dx, b.x);
        T dq = sum_of_products(c.x, d.dx, c.dx, d.x);
        return Dual<T,P>(
            diff_of_products(a, b, c, d),
            dp - dq);
    }
    
    static inline Dual<T,P> sum_of_products(Dual<T,P> a, Dual<T,P> b, Dual<T,P> c, Dual<T,P> d) {
        T dp = sum_of_products(a.x, b.dx, a.dx, b.x);
        T dq = sum_of_products(c.x, d.dx, c.dx, d.x);
        return Dual<T,P>(
            sum_of_products(a, b, c, d),
            dp + dq);
    }
    
    static inline Dual<T,P> multiply_add(Dual<T,P> a, Dual<T,P> b, Dual<T,P> c) {
        // implemented in Dual.h
        return std::fma(a, b, c);
    }
};

}
}

#endif
