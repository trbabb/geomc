/* 
 * File:   Basis.h
 * Author: tbabb
 *
 * Created on December 24, 2014, 3:51 PM
 */

#ifndef BASIS_H
#define	BASIS_H

#include <math.h>
#include <numeric>
#include <boost/math/constants/constants.hpp>
#include <geomc/linalg/Vec.h>
#include <geomc/function/FunctionTypes.h>


namespace geom {

/**
 * @addtogroup function
 * @{
 */
    
// helper function evaluates Legendre polynomials in batch.
// makes use of recurrent definition of legendre polys to save re-execs.
// evaluate a column m and -m of Legendre polynomials and stuff them
// into the appropriate place in the SphericalHarmonics' coeff table
template <typename T, index_t N>
inline void legendre(SphericalHarmonics<T,N> *sh, index_t m, T x) {
    index_t l = sh->bands();
    T pmm0 = 1;
    if (m > 0) {
        // compute double factorial of 2m - 1
        // fold in the (-1)^m
        for (index_t i = 1; i <= m; i++) {
            pmm0 *= -(i * 2 - 1);
        }
        pmm0 *= std::pow(1 - x * x, m / 2.0);
    }
    sh->coeff(m,m) = sh->coeff(m,-m) = pmm0;
    if (m == l - 1) return; // there is no higher band
    T pmm1 = sh->coeff(m+1,m) = sh->coeff(m+1,-m) = x * (2 * m + 1) * pmm0;
    if (m == l - 2) return; // there is no higher band
    for (index_t l_i = m + 2; l_i < l; l_i++) {
        T pmi = (x * (2 * l_i - 1) * pmm1 - (l_i + m - 1) * pmm0) / (l_i - m);
        pmm0 = pmm1;
        pmm1 = sh->coeff(l_i,m) = sh->coeff(l_i,-m) = pmi;
    }
}

/**
 * Evaluate an associated Legendre polynomial.
 * 
 * @code#include <geomc/function/Basis.h>@endcode
 * @param l Band index.
 * @param m Sub-band.
 * @param x Value of `x` at which to evaluate the polynomial.
 */
template <typename T>
inline T legendre(index_t l, index_t m, T x) {
    T pmm0 = 1;
    if (m > 0) {
        // compute double factorial of 2m - 1
        // fold in the (-1)^m
        for (index_t i = 1; i <= m; i++) {
            pmm0 *= -(i * 2 - 1);
        }
        pmm0 *= std::pow(1 - x * x, m / 2.0);
    }
    if (l == m) return pmm0;
    T pmm1 = x * (2 * m + 1) * pmm0;
    if (l == m + 1) return pmm1;
    for (index_t i = m + 2; i <= l; i++) {
        T pmi = (x * (2 * i - 1) * pmm1 - (i + m - 1) * pmm0) / (i - m);
        pmm0 = pmm1;
        pmm1 = pmi;
    }
    return pmm1;
}

/**
 * Evaluate a Legendre polynomal of order `n` at `x`. Equivalent to the associated
 * legendre polynomial with m = 0.
 * 
 * @code#include <geomc/function/Basis.h>@endcode
 * @param n Order of polynomial.
 * @param x Value of `x` at which to evaluate the polynomial.
 */
template <typename T>
inline T legendre(index_t n, T x) {
    if (n <= 0) return 1;
    T p_0 = 1; // P(i - 1, x)
    T p_1 = x; // P(i,     x)
    T p_i = x; // P(i + 1, x)
    for (index_t i = 1; i < n; i++) {
        p_i = ((2 * i + 1) * x * p_1 - i * p_0) / ((T)i + 1);
        p_0 = p_1;
        p_1 = p_i;
    }
    return p_i;
}

/**
 * Evaluate the integral of the Legendre polynomal of order `n` 
 * between -1 and x.
 * 
 * @code#include <geomc/function/Basis.h>@endcode
 * @param n Order of polynomial.
 * @param x Value of `x` at which to evaluate the integral.
 */
template <typename T>
inline T legendre_integral(index_t n, T x) {
    if (n <= 0) return x + 1;
    T p_0 = 1; // P(i - 1, x)
    T p_1 = x; // P(i,     x)
    T p_i = x; // P(i + 1, x)
    for (index_t i = 1; i < n; i++) {
        p_i = ((2 * i + 1) * x * p_1 - i * p_0) / ((T)i + 1);
        p_0 = p_1;
        p_1 = p_i;
    }
    p_i = ((2 * n + 1) * x * p_1 - n * p_0) / ((T)n + 1);
    return (p_i - p_0) / (2 * n + 1);
}

/**
 * Evaluate a Chebyshev polynomial of order `n` at `x`.
 * 
 * @code#include <geomc/function/Basis.h>@endcode
 * @param kind 1 or 2 for the Chebyshev polynomial of the first or second kind, respectively.
 * @param n Order of the polynomial.
 * @param x Value of `x` at which to evaluate.
 */
template <typename T>
inline T chebyshev(index_t kind, index_t n, T x) {
    if (n <= 0) return 1;
    T t_0 = 1;
    T t_1 = kind * x; // first kind := x; second kind := 2x
    T t_i = t_1;
    for (index_t i = 1; i < n; i++) {
        t_i = 2 * x * t_1 - t_0;
        t_0 = t_1;
        t_1 = t_i;
    }
    return t_i;
}


/// @} // group function

} // end namespace geom

#endif	/* BASIS_H */

