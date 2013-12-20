/* 
 * File:   Dual.h
 * Author: tbabb
 *
 * Created on December 18, 2013, 12:50 AM
 */

#ifndef DUAL_H
#define	DUAL_H

#include <geomc/geomc_defs.h>
#include <cmath>
#include <limits>

#ifdef GEOMC_FUNCTION_USE_STREAMS
#include <iostream>
#endif

//todo: what if you have a vector of duals? Does that work, or would the mult operators be hidden?

/* why the comparison operators should NOT be overloaded
 *   should inequality be implemented?
 *     - if so, should epsilon count?
 *     - if so, then < and > and == and != must be in agreement. See below.
 *   should (5, 1) and (5, 0) compare equally?
 *     - if not, then these are not a drop-in replacement for the reals, and errors
 *       are likely to show up because you think they are but aren't forced to think about it.
 *     - if so, then (5,1) == (5,0) is very unexpected. the two have different meaning
 *       and should not be equal.
 */ 

namespace geom {
    
    
/**
 * @ingroup function
 * @{
 */

/**
 * @brief Class implementing the dual numbers, whose arithmetic operations produce a 
 * simultaneous calculation of the first derivative.
 * 
 * Overview
 * ========
 * 
 * Dual numbers implicitly compute and store the first derivative of any arbitrary 
 * function by extending the reals to include a new nonzero element &epsilon; 
 * whose square is zero. Thus duals have the form (a + b&epsilon;); where `a` 
 * may be thought of as the function value, and `b` its first derivative at `a`. 
 * This technique is commonly known as Automatic Differentiation, or AD.
 * 
 * All fundamental operations on duals (including addition, multiplication, division, 
 * and any number of other primitive mathematical functions) implicitly compute and 
 * keep track of the derivative by performing the chain and product rules in-place.
 * Efficiency and accuracy is very good (easily better than either symbolic or 
 * numerical differentiation), adding only a small constant factor to arithmetic
 * operators.
 * 
 * For more on fundamental operations extended to support Duals, see the @ref std
 * "std namespace documentation". 
 * 
 * Use
 * ===
 * 
 * functions: duals to the bottom
 * functions: directional derivatives.
 * creating new fundamental operators.
 * 
 * Discontinuity behavior
 * ======================
 */
template <typename T>
class Dual {
    public:
    
    /// Real component
    T x;
    /// Dual (epsilon) component
    T dx;
    
    /*******************************
     * Constructors                *
     *******************************/
   
    /// Construct a dual with 0 value and 0 derivative.
    Dual<T>():x(0),dx(0) {}
    
    /// Construct a dual with value `x` and 0 derivative.
    Dual<T>(T x):x(x),dx(0) {}
    
    /// Construct a dual with value `x` and derivative `dx`.
    Dual<T>(T x, T dx):x(x),dx(dx) {}
    
    
    /*******************************
     * Operators                   *
     *******************************/
    
    // mult
    
    /// Multiplication.
    friend inline Dual<T> operator*(const Dual<T> &d1, const Dual<T> &d2) {
        return Dual<T>(d1.x * d2.x, d1.x * d2.dx + d1.dx * d2.x);
    }
    
    /// Dual-scalar multiplication.
    template <typename U>
    friend inline Dual<T> operator*(const Dual<T> &d1, U s) {
        return Dual<T>(d1.x * s, d1.dx * s);
    }
    
    /// Scalar-dual multiplication.
    template <typename U>
    friend inline Dual<T> operator*(U s, const Dual<T> &d1) {
        return Dual<T>(s * d1.x, s * d1.dx);
    }
    
    /// Multiply and store.
    Dual<T>& operator*=(const Dual<T> &d) {
        x  = x * d.x;
        dx = x * d.dx + dx * d.x;
        return *this;
    }
    
    /// Multiply and store.
    template <typename U>
    Dual<T>& operator*=(U s) {
        x  =  x * s;
        dx = dx * s;
        return *this;
    }
    
    // div
    
    /// Division.
    friend inline Dual<T> operator/(const Dual<T> &d1, const Dual<T> &d2) {
        return Dual<T>(d1.x / d2.x, (d1.dx * d2.x - d1.x * d2.dx) / (d2.x * d2.x));
    }
    
    /// Dual / scalar division.
    template <typename U>
    friend inline Dual<T> operator/(const Dual<T> &d, U s) {
        return Dual<T>(d.x / s, d.dx / s);
    }
    
    /// Scalar / dual division.
    template <typename U>
    friend inline Dual<T> operator/(U s, const Dual<T> &d) {
        return Dual<T>(s / d.x, -s*d.dx / (d.x * d.x) );
    }
    
    /// Divide and store.
    Dual<T>& operator/=(const Dual<T> &d) {
        x  = x / d.x;
        dx = (dx * d.x - x * d.dx) / (d.x * d.x);
        return *this;
    }
    
    /// Divide and store.
    template <typename U>
    Dual<T>& operator/=(U s) {
        x  =  x / s;
        dx = dx / s;
        return *this;
    }
    
    // add, sub
    
    /// Addition.
    friend inline Dual<T> operator+(const Dual<T> &d1, const Dual<T> &d2) {
        return Dual<T>(d1.x + d2.x, d1.dx + d2.dx);
    }
    
    /// Subtraction.
    friend inline Dual<T> operator-(const Dual<T> &d1, const Dual<T> &d2) {
        return Dual<T>(d1.x - d2.x, d1.dx - d2.dx);
    }
    
    /// Add and store.
    Dual<T>& operator+=(const Dual<T> &d) {
        x  += d.x;
        dx += d.dx;
        return *this;
    }
    
    /// Subtract and store.
    Dual<T>& operator-=(const Dual<T> &d) {
        x  -= d.x;
        dx -= d.dx;
        return *this;
    }
    
    // negation
    
    
    /// Negation.
    inline Dual<T> operator-() {
        return Dual<T>(-x, -dx);
    }
    
    // index
    
    /**
     * @brief Indexing.
     * 
     * Real value is element 0, epsilon value is element 1.
     */
    inline T& operator[](index_t idx) {
        return (&x)[idx];
    }
    
    /**
     * @brief Indexing.
     * 
     * Real value is element 0, epsilon value is element 1.
     */
    inline T operator[](index_t idx) const {
        return (&x)[idx];
    }
    
    // stream
    
#ifdef GEOMC_FUNCTION_USE_STREAMS
    
    /// Stream output.
    friend std::ostream &operator<< (std::ostream &stream, const Dual<T> &d) {
        stream << "(" << d.x << " + " << d.dx << " dx)";
        return stream;
    }
    
#endif

};

/// @} //ingroup function

}; // namespace geom



namespace std {
    
    /**
     * @ingroup function
     * @{
     */
    
    // general form:
    //   f(x), x` * f`(x)
    // or in other words:
    //   f(x), dx * f`(x)
    
    // TODO: atan2
    
    /// Sine function.
    template <typename T>
    inline geom::Dual<T> sin(const geom::Dual<T> &d) {
        return geom::Dual<T>(sin(d.x), d.dx * cos(d.x));
    }
    
    
    /// Cosine function.
    template <typename T>
    inline geom::Dual<T> cos(const geom::Dual<T> &d) {
        return geom::Dual<T>(cos(d.x), -d.dx * sin(d.x));
    }
    
    
    /// Tangent function.
    template <typename T>
    inline geom::Dual<T> tan(const geom::Dual<T> &d) {
        T c = cos(d.x);
        return geom::Dual<T>(tan(d.x), -d.dx / (c*c));
    }
    
    
    /// Arcsin (inverse sine) function.
    template <typename T>
    inline geom::Dual<T> asin(const geom::Dual<T> &d) {
        return geom::Dual<T>(
                asin(d.x),
                d.dx / sqrt(1 - d.x * d.x));
    }
    
    
    /// Arccos (inverse cosine) function.
    template <typename T>
    inline geom::Dual<T> acos(const geom::Dual<T> &d) {
        return geom::Dual<T>(
                acos(d.x),
               -d.dx / sqrt(1 - d.x * d.x));
    }
    
    
    /// Arctan (inverse tangent) function.
    template <typename T>
    inline geom::Dual<T> atan(const geom::Dual<T> &d) {
        return geom::Dual<T>(
                atan(d.x),
                d.dx / (d.x * d.x + 1));
    }
    
    
    /// Exponential (e<sup>x</sup>) function.
    template <typename T>
    inline geom::Dual<T> exp(const geom::Dual<T> &d) {
        T e = exp(d.x);
        return geom::Dual<T>(e, -d.dx * e);
    }
    
    
    /// Returns `base` raised to exponent `xp`.
    template <typename T>
    inline geom::Dual<T> pow(const geom::Dual<T> &base, const geom::Dual<T> &xp) {
        T a_c = pow(base.x, xp.x);
        return geom::Dual<T>(
                a_c, 
                base.dx * xp.x * a_c / xp.x + xp.dx * a_c * log(base.x));
    }
    
    
    /// Returns `base` raised to exponent `xp`.
    template <typename T, typename U>
    inline geom::Dual<T> pow(const geom::Dual<T> &base, U xp) {
        T a_x = pow(base.x, xp);
        return geom::Dual<T>(
                a_x,
                base.dx * xp * a_x / base.x);
    }
    
    
    /// Returns `base` raised to exponent `xp`.
    template <typename T, typename U>
    inline geom::Dual<T> pow(U base, const geom::Dual<T> &xp) {
        T a_x = pow(base, xp.x);
        return geom::Dual<T>(
                a_x,
                xp.dx * a_x * log(base));
    }
    
    
    /// Square root.
    template <typename T>
    inline geom::Dual<T> sqrt(const geom::Dual<T> &d) {
        T sr = sqrt(d.x);
        return geom::Dual<T>(
                sr,
                d.dx / (2*sr));
    }
    
    
    // min, max, ceil, floor
    
    
    template <typename T>
    inline geom::Dual<T> abs(const geom::Dual<T> &d) {
        bool neg = d.x < 0;
        T dx;
        
        if (d.x == 0) {
#if DUAL_DISCONTINUITY_NAN
            dx = numeric_limits<T>::quiet_NaN();
#elif DUAL_DISCONTINUITY_AVERAGE
            dx = 0;
#elif DUAL_DISCONTINUITY_LEFT
            dx = d.dx > 0 ? -d.dx : d.dx;
#elif DUAL_DISCONTINUITY_RIGHT
            dx = d.dx < 0 ? -d.dx : d.dx;
#endif
        } else {
            dx = neg ? -d.dx: d.dx;
        }
        
        return geom::Dual<T>(x,dx);
    }
    
    
    template <typename T>
    inline geom::Dual<T> ceil(const geom::Dual<T> &d) {
        T x = ceil(d.x);
#if DUAL_DISCONTINUITY_NAN
        T dx = (x == d.x) ? numeric_limits<T>::quiet_NaN() : 0;
#else
        T dx = 0;
#endif
        return geom::Dual<T>(x,dx);
    }
    
    
    template <typename T>
    inline geom::Dual<T> floor(const geom::Dual<T> &d) {
        T x = floor(d.x);
#if DUAL_DISCONTINUITY_NAN
        T dx = (x == d.x) ? numeric_limits<T>::quiet_NaN() : 0;
#else
        T dx = 0;
#endif
        return geom::Dual<T>(x,dx);
    }
    
    
    template <typename T>
    inline geom::Dual<T> min(const geom::Dual<T> &d1, const geom::Dual<T> &d2) {
        bool lt = d1.x <= d2.x;
        T x = lt ? d1.x : d2.x;
        T dx;
        if (d1.x == d2.x) {
#if DUAL_DISCONTINUITY_NAN
            dx = numeric_limits<T>::quiet_NaN();
#elif DUAL_DISCONTINUITY_AVERAGE
            dx = (d1.dx + d2.dx) / 2;
#elif DUAL_DISCONTINUITY_LEFT
            dx = (d1.dx < d2.dx) ? d1.dx : d2.dx;
#elif DUAL_DISCONTINUITY_RIGHT
            dx = (d1.dx > d2.dx) ? d1.dx : d2.dx;
#endif
        } else {
            dx = lt ? d1.dx : d2.dx;
        }
        return geom::Dual<T>(x,dx);
    }
    
        
    template <typename T>
    inline geom::Dual<T> max(const geom::Dual<T> &d1, const geom::Dual<T> &d2) {
        bool gt = d1.x >= d2.x;
        T x = gt ? d1.x : d2.x;
        T dx;
        if (d1.x == d2.x) {
#if DUAL_DISCONTINUITY_NAN
            dx = numeric_limits<T>::quiet_NaN();
#elif DUAL_DISCONTINUITY_AVERAGE
            dx = (d1.dx + d2.dx) / 2;
#elif DUAL_DISCONTINUITY_LEFT
            dx = (d1.dx > d2.dx) ? d1.dx : d2.dx;
#elif DUAL_DISCONTINUITY_RIGHT
            dx = (d1.dx < d2.dx) ? d1.dx : d2.dx;
#endif
        } else {
            dx = gt ? d1.dx : d2.dx;
        }
        return geom::Dual<T>(x,dx);
    }

    /// @} // ingroup function
    
};

#endif	/* DUAL_H */

