/* 
 * File:   Dual.h
 * Author: tbabb
 *
 * Created on December 18, 2013, 12:50 AM
 */

#ifndef DUAL_H
#define	DUAL_H

//TODO: FIXME
#include <geomc/geomc_defs.h>
#include <cmath>
#include <limits>

//todo: what if you have a vector of duals? Does that work, or would the mult operators be hidden?

namespace geom {

template <typename T>
class Dual {
    public:
    
    /// Real component
    T a;
    /// Dual (epsilon) component
    T b;
    
    /*******************************
     * Constructors                *
     *******************************/
   
    Dual<T>():a(0),b(0) {}
    
    Dual<T>(T a):a(a),b(0) {}
    
    Dual<T>(T a, T b):a(a),b(b) {}
    
    
    /*******************************
     * Operators                   *
     *******************************/
    
    // mult
    
    friend inline Dual<T> operator*(const Dual<T> &d1, const Dual<T> &d2) {
        return Dual<T>(d1.a * d2.a, d1.a * d2.b + d1.b * d2.a);
    }
    
    template <typename U>
    Dual<T> operator*(const Dual<T> &d1, U s) {
        return Dual<T>(d1.a * s, d1.b * s);
    }
    
    template <typename U>
    Dual<T> operator*(U s, const Dual<T> &d1) {
        return Dual<T>(s * d1.a, s * d1.b);
    }
    
    Dual<T>& operator*=(const Dual<T> &d) {
        a = a * d.a;
        b = a * d.b + b * d.a;
        return *this;
    }
    
    template <typename U>
    Dual<T>& operator*=(U s) {
        a = a * s;
        b = b * s;
        return *this;
    }
    
    // div
    
    friend inline Dual<T> operator/(const Dual<T> &d1, const Dual<T> &d2) {
        return Dual<T>(d1.a / d2.a, (d1.b * d2.a - d1.a * d2.b) / (d2.a * d2.a));
    }
    
    template <typename U>
    friend inline Dual<T> operator/(const Dual<T> &d, U s) {
        return Dual<T>(d.a / s, d.b / s);
    }
    
    template <typename U>
    friend inline Dual<T> operator/(U s, const Dual<T> &d) {
        return Dual<T>(s / d.a, -s*d.b / (d.a * d.a) );
    }
    
    Dual<T>& operator/=(const Dual<T> &d) {
        a = a / d.a;
        b = (b * d.a - a * d.b) / (d.a * d.a);
        return *this;
    }
    
    template <typename U>
    Dual<T>& operator/=(U s) {
        a = a / s;
        b = b / s;
        return *this;
    }
    
    // add, sub
    
    friend inline Dual<T> operator+(const Dual<T> &d1, const Dual<T> &d2) {
        return Dual<T>(d1.a + d2.a, d1.b + d2.b);
    }
    
    friend inline Dual<T> operator-(const Dual<T> &d1, const Dual<T> &d2) {
        return Dual<T>(d1.a - d2.a, d1.b - d2.b);
    }
    
    Dual<T>& operator+=(const Dual<T> &d) {
        a += d.a;
        b += d.b;
        return *this;
    }
    
    Dual<T>& operator-=(const Dual<T> &d) {
        a -= d.a;
        b -= d.b;
        return *this;
    }
    
    // negation
    
    inline Dual<T> operator-() {
        return Dual<T>(-a, -b);
    }
    
    // index
    
    inline T& operator[](index_t idx) {
        return (&a)[idx];
    }
    
    inline T operator[](index_t idx) const {
        return (&a)[idx];
    }
};

}; // namespace geom



namespace std {
    
    // general form:
    //   f(x), x` * f`(x)
    // or in other words:
    //   f(a), b * f`(a)
    
    // TODO: atan2
    
    template <typename T>
    inline geom::Dual<T> sin(const geom::Dual<T> &d) {
        return geom::Dual<T>(sin(d.a), d.b * cos(d.a));
    }
    
    
    template <typename T>
    inline geom::Dual<T> cos(const geom::Dual<T> &d) {
        return geom::Dual<T>(cos(d.a), -d.b * sin(d.a));
    }
    
    
    template <typename T>
    inline geom::Dual<T> tan(const geom::Dual<T> &d) {
        T c = cos(d.a);
        return geom::Dual<T>(tan(d.a), -d.b / (c*c));
    }
    
    
    template <typename T>
    inline geom::Dual<T> asin(const geom::Dual<T> &d) {
        return geom::Dual<T>(
                asin(d.a),
                d.b / sqrt(1 - d.a * d.a));
    }
    
    template <typename T>
    inline geom::Dual<T> acos(const geom::Dual<T> &d) {
        return geom::Dual<T>(
                acos(d.a),
               -d.b / sqrt(1 - d.a * d.a));
    }
    
    
    template <typename T>
    inline geom::Dual<T> atan(const geom::Dual<T> &d) {
        return geom::Dual<T>(
                atan(d.a),
                d.b / (d.a * d.a + 1));
    }
    
    
    template <typename T>
    inline geom::Dual<T> exp(const geom::Dual<T> &d) {
        T e = exp(d.a);
        return geom::Dual<T>(e, -d.b * e);
    }
    
    
    template <typename T>
    inline geom::Dual<T> pow(const geom::Dual<T> &base, const geom::Dual<T> &xp) {
        T a_c = pow(base.a, xp.a);
        return geom::Dual<T>(
                a_c, 
                base.b * xp.a * a_c / xp.a + xp.b * a_c * log(base.a));
    }
    
    
    template <typename T, typename U>
    inline geom::Dual<T> pow(const geom::Dual<T> &base, U xp) {
        T a_x = pow(base.a, xp);
        return geom::Dual<T>(
                a_x,
                base.b * xp * a_x / base.a);
    }
    
    
    template <typename T, typename U>
    inline geom::Dual<T> pow(U base, const geom::Dual<T> &xp) {
        T a_x = pow(base, xp.a);
        return geom::Dual<T>(
                a_x,
                xp.b * a_x * log(base));
    }
    
    
    template <typename T>
    inline geom::Dual<T> sqrt(const geom::Dual<T> &d) {
        T sr = sqrt(d.a);
        return geom::Dual<T>(
                sr,
                d.b / (2*sr));
    }
    
    
    // min, max, ceil, floor
    
    
    template <typename T>
    inline geom::Dual<T> abs(const geom::Dual<T> &d) {
#if DUAL_DISCONTINUITY_LEFT
        bool pos = d.a > 0;
#else
        bool pos = d.a >= 0;
#endif
        T a = pos ? d.a : -d.a;
        T b;
#if !(DUAL_DISCONTINUITY_LEFT || DUAL_DISCONTINUITY_RIGHT)
        if (a == 0) {
#if DUAL_DISCONTINUITY_NAN
            a = numeric_limits<T>::quiet_NaN();
#elif DUAL_DISCONTINUITY_AVERAGE
            a = 0;
#endif
        } else
#endif
        b = pos ? d.b : -d.b;
        
        return geom::Dual<T>(a,b);
    }
    
    
    template <typename T>
    inline geom::Dual<T> ceil(const geom::Dual<T> &d) {
        T a = ceil(d.a);
#if DUAL_DISCONTINUITY_NAN
        T b = (a == d.a) ? numeric_limits<T>::quiet_NaN() : 0;
#else
        T b = 0;
#endif
        return geom::Dual<T>(a,b);
    }
    
    
    template <typename T>
    inline geom::Dual<T> floor(const geom::Dual<T> &d) {
        T a = floor(d.a);
#if DUAL_DISCONTINUITY_NAN
        T b = (a == d.a) ? numeric_limits<T>::quiet_NaN() : 0;
#else
        T b = 0;
#endif
        return geom::Dual<T>(a,b);
    }
    
    
    template <typename T>
    inline geom::Dual<T> min(const geom::Dual<T> &d1, const geom::Dual<T> &d2) {
#if DUAL_DISCONTINUITY_LEFT
        bool lt = d1.a <= d2.a;
#else
        bool lt = d1.a < d2.a;
#endif
        T a = lt ? d1.a : d2.a;
        T b;
#if DUAL_DISCONTINUITY_NAN || DUAL_DISCONTINUITY_AVERAGE
        if (d1.a == d2.a) {
#if DUAL_DISCONTINUITY_NAN
            b = numeric_limits<T>::quiet_NaN();
#else
            b = (d1.b + d2.b) / 2
#endif
        } else 
#endif
        b = lt ? d1.b : d2.b;
        return geom::Dual<T>(a,b);
    }
    
        
    template <typename T>
    inline geom::Dual<T> max(const geom::Dual<T> &d1, const geom::Dual<T> &d2) {
#if DUAL_DISCONTINUITY_LEFT
        bool gt = d1.a >= d2.a;
#else
        bool gt = d1.a > d2.a;
#endif
        T a = gt ? d1.a : d2.a;
        T b;
#if DUAL_DISCONTINUITY_NAN || DUAL_DISCONTINUITY_AVERAGE
        if (d1.a == d2.a) {
#if DUAL_DISCONTINUITY_NAN
            b = numeric_limits<T>::quiet_NaN();
#else
            b = (d1.b + d2.b) / 2
#endif
        } else 
#endif
        b = gt ? d1.b : d2.b;
        return geom::Dual<T>(a,b);
    }
    
};

#endif	/* DUAL_H */

