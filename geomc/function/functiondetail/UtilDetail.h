#pragma once

#include <cmath>
#include <type_traits>
#include <complex>

#include <geomc/function/functiondetail/DualImpl.h>

namespace geom {

///////////////
// fwd decls //
///////////////

template <typename T>
inline T multiply_add(T a, T b, T c);

template <typename T>
inline T diff_of_products(T a, T b, T c, T d);

template <typename T>
inline T sum_of_products(T a, T b, T c, T d);


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
struct _ImplFMADot<T, typename std::enable_if_t<std::is_floating_point_v<T>>> {
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
    static inline Dual<T,P> diff_of_products(
            Dual<T,P> a,
            Dual<T,P> b,
            Dual<T,P> c,
            Dual<T,P> d)
    {
        T dp = geom::sum_of_products<T>(a.x, b.dx, a.dx, b.x);
        T dq = geom::sum_of_products<T>(c.x, d.dx, c.dx, d.x);
        return Dual<T,P>(
            geom::diff_of_products<T>(a, b, c, d),
            dp - dq);
    }
    
    static inline Dual<T,P> sum_of_products(
            Dual<T,P> a, 
            Dual<T,P> b, 
            Dual<T,P> c, 
            Dual<T,P> d)
    {
        T dp = geom::sum_of_products<T>(a.x, b.dx, a.dx, b.x);
        T dq = geom::sum_of_products<T>(c.x, d.dx, c.dx, d.x);
        return Dual<T,P>(
            geom::sum_of_products<T>(a, b, c, d),
            dp + dq);
    }
    
    static inline Dual<T,P> multiply_add(Dual<T,P> a, Dual<T,P> b, Dual<T,P> c) {
        return Dual<T,P>(
            geom::multiply_add<T>(a.x, b.x, c.x),
            geom::sum_of_products<T>(a.x, b.dx, a.dx, b.x) + c.dx
        );
    }
};


template <typename T>
struct _ImplFMADot<std::complex<T>, void> {
    static inline std::complex<T> diff_of_products(
            std::complex<T> a,
            std::complex<T> b,
            std::complex<T> c,
            std::complex<T> d)
    {
        // (a.r * b.r - a.i * b.i) + (a.r * b.i + a.i * b.r)i -
        // (c.r * d.r - c.i * d.i) + (c.r * d.i + c.i * d.r)i
        auto ab_r = geom::diff_of_products<T>(a.real, b.real, a.imag, b.imag);
        auto cd_r = geom::diff_of_products<T>(c.real, d.real, c.imag, d.imag);
        auto ab_i =  geom::sum_of_products<T>(a.real, b.imag, a.imag, b.real);
        auto cd_i =  geom::sum_of_products<T>(c.real, d.imag, c.imag, d.real);
        return {ab_r - cd_r, ab_i - cd_i};
    }
    
    static inline std::complex<T> sum_of_products(
            std::complex<T> a,
            std::complex<T> b,
            std::complex<T> c,
            std::complex<T> d)
    {
        // (a.r * b.r - a.i * b.i) + (a.r * b.i + a.i * b.r)i +
        // (c.r * d.r - c.i * d.i) + (c.r * d.i + c.i * d.r)i
        auto ab_r = geom::diff_of_products<T>(a.real, b.real, a.imag, b.imag);
        auto cd_r = geom::diff_of_products<T>(c.real, d.real, c.imag, d.imag);
        auto ab_i =  geom::sum_of_products<T>(a.real, b.imag, a.imag, b.real);
        auto cd_i =  geom::sum_of_products<T>(c.real, d.imag, c.imag, d.real);
        return {ab_r + cd_r, ab_i + cd_i};
    }
    
    static inline std::complex<T> multiply_add(
            std::complex<T> a, 
            std::complex<T> b, 
            std::complex<T> c)
    {
        // todo: better...?
        // (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
        // auto ab_r = diff_of_products(a.real, b.real, a.imag, b.imag);
        // auto ab_i =  sum_of_products(a.real, b.imag, a.imag, b.real);
        // return {ab_r + c.real, ab_i + c.i};
        
        // or maybe try to kahanify this yourself...?
        
        // (a + bi) * (c + di) + (k_r + k_i i)
        // real: (ac - bd + k_r)
        auto bd_k = geom::multiply_add<T>(-a.imag, b.imag, c.real);
        auto real = geom::multiply_add<T>( a.real, b.real, bd_k);
        // imag: (ad + bc + k_i)
        auto bc_k = geom::multiply_add<T>(a.imag, b.real, c.imag);
        auto imag = geom::multiply_add<T>(a.real, b.imag, bc_k);
        return {real, imag};
    }
};

} // detail
} // geom
