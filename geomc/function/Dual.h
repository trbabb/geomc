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

namespace geom {

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
   
    Dual<T>():x(0),dx(0) {}
    
    Dual<T>(T x):x(x),dx(0) {}
    
    Dual<T>(T x, T dx):x(x),dx(dx) {}
    
    
    /*******************************
     * Operators                   *
     *******************************/
    
    // mult
    
    friend inline Dual<T> operator*(const Dual<T> &d1, const Dual<T> &d2) {
        return Dual<T>(d1.x * d2.x, d1.x * d2.dx + d1.dx * d2.x);
    }
    
    template <typename U>
    friend inline Dual<T> operator*(const Dual<T> &d1, U s) {
        return Dual<T>(d1.x * s, d1.dx * s);
    }
    
    template <typename U>
    friend inline Dual<T> operator*(U s, const Dual<T> &d1) {
        return Dual<T>(s * d1.x, s * d1.dx);
    }
    
    Dual<T>& operator*=(const Dual<T> &d) {
        x  = x * d.x;
        dx = x * d.dx + dx * d.x;
        return *this;
    }
    
    template <typename U>
    Dual<T>& operator*=(U s) {
        x  =  x * s;
        dx = dx * s;
        return *this;
    }
    
    // div
    
    friend inline Dual<T> operator/(const Dual<T> &d1, const Dual<T> &d2) {
        return Dual<T>(d1.x / d2.x, (d1.dx * d2.x - d1.x * d2.dx) / (d2.x * d2.x));
    }
    
    template <typename U>
    friend inline Dual<T> operator/(const Dual<T> &d, U s) {
        return Dual<T>(d.x / s, d.dx / s);
    }
    
    template <typename U>
    friend inline Dual<T> operator/(U s, const Dual<T> &d) {
        return Dual<T>(s / d.x, -s*d.dx / (d.x * d.x) );
    }
    
    Dual<T>& operator/=(const Dual<T> &d) {
        x  = x / d.x;
        dx = (dx * d.x - x * d.dx) / (d.x * d.x);
        return *this;
    }
    
    template <typename U>
    Dual<T>& operator/=(U s) {
        x  =  x / s;
        dx = dx / s;
        return *this;
    }
    
    // add, sub
    
    friend inline Dual<T> operator+(const Dual<T> &d1, const Dual<T> &d2) {
        return Dual<T>(d1.x + d2.x, d1.dx + d2.dx);
    }
    
    friend inline Dual<T> operator-(const Dual<T> &d1, const Dual<T> &d2) {
        return Dual<T>(d1.x - d2.x, d1.dx - d2.dx);
    }
    
    Dual<T>& operator+=(const Dual<T> &d) {
        x  += d.x;
        dx += d.dx;
        return *this;
    }
    
    Dual<T>& operator-=(const Dual<T> &d) {
        x  -= d.x;
        dx -= d.dx;
        return *this;
    }
    
    // negation
    
    inline Dual<T> operator-() {
        return Dual<T>(-x, -dx);
    }
    
    // index
    
    inline T& operator[](index_t idx) {
        return (&x)[idx];
    }
    
    inline T operator[](index_t idx) const {
        return (&x)[idx];
    }
    
    // stream
    
#ifdef GEOMC_FUNCTION_USE_STREAMS
    
    friend std::ostream &operator<< (std::ostream &stream, const Dual<T> &d) {
        stream << "(" << d.x << " + " << d.dx << " dx)";
        return stream;
    }
    
#endif

};

}; // namespace geom



namespace std {
    
    // general form:
    //   f(x), x` * f`(x)
    // or in other words:
    //   f(x), dx * f`(x)
    
    // TODO: atan2
    
    template <typename T>
    inline geom::Dual<T> sin(const geom::Dual<T> &d) {
        return geom::Dual<T>(sin(d.x), d.dx * cos(d.x));
    }
    
    
    template <typename T>
    inline geom::Dual<T> cos(const geom::Dual<T> &d) {
        return geom::Dual<T>(cos(d.x), -d.dx * sin(d.x));
    }
    
    
    template <typename T>
    inline geom::Dual<T> tan(const geom::Dual<T> &d) {
        T c = cos(d.x);
        return geom::Dual<T>(tan(d.x), -d.dx / (c*c));
    }
    
    
    template <typename T>
    inline geom::Dual<T> asin(const geom::Dual<T> &d) {
        return geom::Dual<T>(
                asin(d.x),
                d.dx / sqrt(1 - d.x * d.x));
    }
    
    template <typename T>
    inline geom::Dual<T> acos(const geom::Dual<T> &d) {
        return geom::Dual<T>(
                acos(d.x),
               -d.dx / sqrt(1 - d.x * d.x));
    }
    
    
    template <typename T>
    inline geom::Dual<T> atan(const geom::Dual<T> &d) {
        return geom::Dual<T>(
                atan(d.x),
                d.dx / (d.x * d.x + 1));
    }
    
    
    template <typename T>
    inline geom::Dual<T> exp(const geom::Dual<T> &d) {
        T e = exp(d.x);
        return geom::Dual<T>(e, -d.dx * e);
    }
    
    
    template <typename T>
    inline geom::Dual<T> pow(const geom::Dual<T> &base, const geom::Dual<T> &xp) {
        T a_c = pow(base.x, xp.x);
        return geom::Dual<T>(
                a_c, 
                base.dx * xp.x * a_c / xp.x + xp.dx * a_c * log(base.x));
    }
    
    
    template <typename T, typename U>
    inline geom::Dual<T> pow(const geom::Dual<T> &base, U xp) {
        T a_x = pow(base.x, xp);
        return geom::Dual<T>(
                a_x,
                base.dx * xp * a_x / base.x);
    }
    
    
    template <typename T, typename U>
    inline geom::Dual<T> pow(U base, const geom::Dual<T> &xp) {
        T a_x = pow(base, xp.x);
        return geom::Dual<T>(
                a_x,
                xp.dx * a_x * log(base));
    }
    
    
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
#if DUAL_DISCONTINUITY_LEFT
        bool pos = d.x > 0;
#else
        bool pos = d.x >= 0;
#endif
        T x = pos ? d.x : -d.x;
        T dx;
#if !(DUAL_DISCONTINUITY_LEFT || DUAL_DISCONTINUITY_RIGHT)
        if (x == 0) {
#if DUAL_DISCONTINUITY_NAN
            x = numeric_limits<T>::quiet_NaN();
#elif DUAL_DISCONTINUITY_AVERAGE
            x = 0;
#endif
        } else
#endif
        dx = pos ? d.dx : -d.dx;
        
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
#if DUAL_DISCONTINUITY_LEFT
        bool lt = d1.x <= d2.x;
#else
        bool lt = d1.x < d2.x;
#endif
        T x = lt ? d1.x : d2.x;
        T dx;
#if DUAL_DISCONTINUITY_NAN || DUAL_DISCONTINUITY_AVERAGE
        if (d1.x == d2.x) {
#if DUAL_DISCONTINUITY_NAN
            dx = numeric_limits<T>::quiet_NaN();
#else
            dx = (d1.dx + d2.dx) / 2
#endif
        } else 
#endif
        dx = lt ? d1.dx : d2.dx;
        return geom::Dual<T>(x,dx);
    }
    
        
    template <typename T>
    inline geom::Dual<T> max(const geom::Dual<T> &d1, const geom::Dual<T> &d2) {
#if DUAL_DISCONTINUITY_LEFT
        bool gt = d1.x >= d2.x;
#else
        bool gt = d1.x > d2.x;
#endif
        T x = gt ? d1.x : d2.x;
        T dx;
#if DUAL_DISCONTINUITY_NAN || DUAL_DISCONTINUITY_AVERAGE
        if (d1.x == d2.x) {
#if DUAL_DISCONTINUITY_NAN
            dx = numeric_limits<T>::quiet_NaN();
#else
            dx = (d1.dx + d2.dx) / 2
#endif
        } else 
#endif
        dx = gt ? d1.dx : d2.dx;
        return geom::Dual<T>(x,dx);
    }
    
};

#endif	/* DUAL_H */

