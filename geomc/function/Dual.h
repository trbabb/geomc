/* 
 * File:   Dual.h
 * Author: tbabb
 *
 * Created on December 18, 2013, 12:50 AM
 */

// Notice: This #define is referred to in other places, namely Random.h.
//         This is part of a mechanism to bring in code which is common
//         to both `random` and `function` iff both libraries are in use.
#ifndef GEOMC_DUAL_H
#define	GEOMC_DUAL_H

#include <type_traits>

#include <geomc/geomc_defs.h>
#include <geomc/Hash.h>
#include <cmath>
#include <limits>

#ifdef GEOMC_FUNCTION_USE_STREAMS
#include <iostream>
#endif

// todo: what if you have a vector of duals? 
//   Does that work, or would the mult operators be hidden?
// todo: extend numeric limits?
// todo: c++11 needed for conversions to Vec<index_t,N> from Vec<Dual,N>
// todo: c++11 needed for vectors of non-POD classes. 
// todo: don't use "scalar", use convertible_to

namespace geom {

/**
 * @brief Policy for handling discontinuous derivatives
 *
 * @ingroup function
 */
enum class DiscontinuityPolicy {
    /// Return `NaN` at discontinuities
    NaN,
    /// At discontinuities, return the average of the two boundary values.
    Average,
    /// At discontinuities, return the value when approaching from the left
    Left,
    /// At discontinuities, return the value when approaching from the right
    Right,
    /**
     * @brief At discontinuities, return +/- infinity, according to whether the
     * discontinuity is an increasing or decreasing jump.
     */
    Inf
};
    
namespace detail {

template <typename T, DiscontinuityPolicy P>
struct Discontinuity {};

template <typename T>
struct Discontinuity<T, DiscontinuityPolicy::NaN> {
    static inline T resolve(T a, T b) {
        return std::numeric_limits<T>::quiet_NaN();
    }
};

template <typename T>
struct Discontinuity<T, DiscontinuityPolicy::Left> {
    static inline T resolve(T a, T b) {
        return a;
    }
};

template <typename T>
struct Discontinuity<T, DiscontinuityPolicy::Right> {
    static inline T resolve(T a, T b) {
        return b;
    }
};

template <typename T>
struct Discontinuity<T, DiscontinuityPolicy::Average> {
    static inline T resolve(T a, T b) {
        return (a + b) / (T)2;
    }
};

template <typename T>
struct Discontinuity<T, DiscontinuityPolicy::Inf> {
    static inline T resolve(T a, T b) {
        const T inf = std::numeric_limits<T>::infinity();
        return (b == a) ? a : (b > a ? inf : -inf);
    }
};

} // namespace detail


/**
 * @ingroup function
 * @{
 */

/**
 * @brief Class implementing the dual numbers, whose arithmetic operations perform a 
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
 * Fundamental operations on duals (including addition, multiplication, division, 
 * and any number of other primitive mathematical functions) implicitly compute and 
 * keep track of the derivative by performing the chain and product rules in-place.
 * Efficiency and accuracy is very good (easily better than either symbolic or 
 * numerical differentiation), adding only a constant factor to arithmetic
 * operations.
 * 
 * In general, Duals should behave exactly as bare numbers in code.
 *
 * It is safe to convert a constant value to a Dual (and this will be done implicitly
 * by setting the epsilon component to zero). It is unsafe to convert a Dual number
 * to a bare value, as this forgets the derivative information. Therefore if the
 * Dual is to be used as a constant, it must be cast explicitly, or the primal
 * component must be extracted.
 * 
 * Use
 * ===
 * 
 * The value of an arbitrary function may be differentiated with respect
 * to its input if and only if those inputs, and all the calculations that depend on
 * those inputs, are performed entirely on Dual numbers. In other words, a function
 * is differentiable if it uses Duals "all the way to the bottom" of its calculation.
 * If ever the primal is extracted and used in some sub-calculation apart from its
 * epsilon component, then the resultant derivative information will be invalid.
 * 
 * The epsilon components of a function's input or output together represent the
 * components of a directional derivative.
 * 
 * Equality
 * ========
 *
 * Duals are compared only using their primal value, so that they behave just as
 * bare numbers do. This means that two `Dual`s with different epsilons may compare
 * equal! If you need exact binary equality, then explicitly check that both `x`s
 * and `dx`s are equal. 
 *
 * @tparam T The numerical type of the number.
 * @tparam DP The `DiscontinuityPolicy` for handling discontinuous derivatives.
 */
template <typename T, DiscontinuityPolicy DP=DiscontinuityPolicy::Right>
class Dual {
    public:
    typedef detail::Discontinuity<T, DP> discontinuity_policy;
    
    /// Real (primal) component
    T x;
    /// Dual (epsilon) component
    T dx;
    
    /*******************************
     * Constructors                *
     *******************************/
   
    /// Construct a dual with 0 value and 0 derivative.
    Dual():x(0),dx(0) {}
    
    /// Construct a dual with value `x` and 0 derivative.
    Dual(T x):x(x),dx(0) {}
    
    /// Construct a dual with value `x` and derivative `dx`.
    Dual(T x, T dx):x(x),dx(dx) {}
    
    template <DiscontinuityPolicy DP1>
    Dual(const Dual<T,DP1>& d):x(d.dx), dx(d.dx) {}
    
    
    /*******************************
     * Operators                   *
     *******************************/
    
    // mult
    
    /// Multiplication.
    friend inline Dual<T> operator*(const Dual<T> &d1, const Dual<T> &d2) {
        return Dual<T>(
            d1.x * d2.x, 
            sum_of_products(d1.x, d2.dx, d1.dx, d2.x));
    }
    
    /// Dual-scalar multiplication.
    template <typename U>
    friend inline typename std::enable_if<std::is_scalar<U>::value, Dual<T>>::type 
    operator*(const Dual<T> &d1, U s) {
        return Dual<T>(d1.x * s, d1.dx * s);
    }
    
    /// Scalar-dual multiplication.
    template <typename U>
    friend inline typename std::enable_if<std::is_scalar<U>::value, Dual<T>>::type 
    operator*(U s, const Dual<T> &d1) {
        return Dual<T>(s * d1.x, s * d1.dx);
    }
    
    /// Multiply and assign.
    Dual<T>& operator*=(const Dual<T> &d) {
        x  = x * d.x;
        dx = x * d.dx + dx * d.x;
        return *this;
    }
    
    /// Multiply and assign.
    template <typename U>
    inline typename std::enable_if<std::is_scalar<U>::value, Dual<T>&>::type  
    operator*=(U s) {
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
    friend inline typename std::enable_if<std::is_scalar<U>::value, Dual<T> >::type 
    operator/(const Dual<T> &d, U s) {
        return Dual<T>(d.x / s, d.dx / s);
    }
    
    /// Scalar / dual division.
    template <typename U>
    friend inline typename std::enable_if<std::is_scalar<U>::value, Dual<T> >::type 
    operator/(U s, const Dual<T> &d) {
        return Dual<T>(s / d.x, -s*d.dx / (d.x * d.x) );
    }
    
    /// Divide and assign.
    Dual<T>& operator/=(const Dual<T> &d) {
        x  = x / d.x;
        dx = (dx * d.x - x * d.dx) / (d.x * d.x);
        return *this;
    }
    
    /// Divide and assign.
    template <typename U>
    inline typename std::enable_if<std::is_scalar<U>::value, Dual<T>&>::type  
    operator/=(U s) {
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
    
    /// Add and assign.
    Dual<T>& operator+=(const Dual<T> &d) {
        x  += d.x;
        dx += d.dx;
        return *this;
    }
    
    /// Subtract and assign.
    Dual<T>& operator-=(const Dual<T> &d) {
        x  -= d.x;
        dx -= d.dx;
        return *this;
    }
    
    // negation
    
    
    /// Negation.
    inline Dual<T> operator-() const {
        return Dual<T>(-x, -dx);
    }
    
    // index
    
    /**
     * @brief Indexing.
     * 
     * Real value is element 0, epsilon value is element 1.
     */
    inline T& operator[](index_t idx) { return (&x)[idx]; }
    
    /**
     * @brief Indexing.
     * 
     * Real value is element 0, epsilon value is element 1.
     */
    inline T operator[](index_t idx) const { return (&x)[idx]; }
    
    // comparison
    
    /// Equality of primal component
    inline bool operator==(const Dual<T> &d) { return x == d.x; }
    /// Inequality of primal component
    inline bool operator!=(const Dual<T> &d) { return x != d.x; }
    /// Greater-or-equal-to comparison of primal component
    inline bool operator>=(const Dual<T> &d) { return x >= d.x; }
    /// Less-or-equal-to comparison of primal component
    inline bool operator<=(const Dual<T> &d) { return x <= d.x; }
    /// Greater-than comparison of primal component
    inline bool operator> (const Dual<T> &d) { return x >  d.x; }
    /// Less-than comparison of primal component
    inline bool operator< (const Dual<T> &d) { return x <  d.x; }
    
    // conversion
    
    template <typename U>
    inline operator Dual<U>() const {
        return Dual<U>((U)x, (U)dx);
    }
    
#if __cplusplus >= 201103L
    
    // older versions of c++ do not support explicit conversion operators.
    // we don't want this to be implicit because the compiler could readily kill 
    // the epsilon component without us knowing.
    template <typename U>
    explicit inline operator U() const { return (U)x; }
    
#endif
    
    // stream
    
#ifdef GEOMC_FUNCTION_USE_STREAMS
    
    /// Stream output.
    friend std::ostream &operator<<(std::ostream& stream, const Dual<T>& d) {
        stream << "(" << d.x;
        if (d.dx < 0) {
            stream << " - " << std::abs(d.dx);
        } else {
            stream << " + " << d.dx;
        }
        stream << " dx)";
        return stream;
    }
    
#endif

}; // class Dual

/// @} //ingroup function


// fwd decl these guys:

template <typename T> 
inline T diff_of_products(T a, T b, T c, T d);

template <typename T>
inline T sum_of_products(T a, T b, T c, T d);

template <typename T>
inline T multiply_add(T a, T b, T c);


} // namespace geom



namespace std {
    
// general form:
//   f(x), x` * f`(x)
// or in other words:
//   f(x), dx * f`(x)

// TODO: atan2

template <class T, geom::DiscontinuityPolicy P>
struct numeric_limits<geom::Dual<T,P>> : std::numeric_limits<T> {};

/**
 * @brief Sine function.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> sin(const geom::Dual<T,P> &d) {
    return geom::Dual<T,P>(sin(d.x), d.dx * cos(d.x));
}


/**
 * @brief Cosine function.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> cos(const geom::Dual<T,P> &d) {
    return geom::Dual<T,P>(cos(d.x), -d.dx * sin(d.x));
}


/**
 * @brief Tangent function.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> tan(const geom::Dual<T,P> &d) {
    T c = cos(d.x);
    return geom::Dual<T,P>(tan(d.x), -d.dx / (c*c));
}


/**
 * @brief Arcsin (inverse sine) function.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> asin(const geom::Dual<T,P> &d) {
    return geom::Dual<T,P>(
            asin(d.x),
            d.dx / sqrt(1 - d.x * d.x));
}


/**
 * @brief Arccos (inverse cosine) function.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> acos(const geom::Dual<T,P> &d) {
    return geom::Dual<T,P>(
            acos(d.x),
            -d.dx / sqrt(1 - d.x * d.x));
}


/**
 * @brief Arctan (inverse tangent) function.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> atan(const geom::Dual<T,P> &d) {
    return geom::Dual<T,P>(
            atan(d.x),
            d.dx / (d.x * d.x + 1));
}


/** 
 * @brief Exponential (e<sup>x</sup>) function.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> exp(const geom::Dual<T,P> &d) {
    T e = exp(d.x);
    return geom::Dual<T,P>(e, -d.dx * e);
}


/** 
 * @brief Returns `base` raised to exponent `xp`.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> pow(const geom::Dual<T,P> &base, const geom::Dual<T,P> &xp) {
    // here we use the chain rule for partial derivatives:
    //   f(x,y) = x ^ y
    //   x = g(t); y = h(t)
    // then:
    //   df / dt = (df / dx) * (dx / dt) + (df / dy) * (dy / dt)
    //           = (df / dx) * g`(t)     + (df / dy) * h`(t)
    // and:
    //   df / dx = d/dx (x ^ y) = y * x ^ (y - 1)
    //   df / dy = d/dy (x ^ y) = log(x) * (x ^ y)
    T a_c = pow(base.x, xp.x);
    return geom::Dual<T,P>(
            a_c, 
            base.dx * xp.x * a_c / xp.x +
            xp.dx * a_c * log(base.x));
    
}


/** 
 * @brief Returns `base` raised to exponent `xp`.
 * @related geom::Dual
 */
template <typename T, typename U>
inline geom::Dual<T> pow(const geom::Dual<T> &base, U xp) {
    T a_x = pow(base.x, xp);
    return geom::Dual<T>(
            a_x,
            base.dx * xp * a_x / base.x);
}


/** 
 * @brief Returns `base` raised to exponent `xp`.
 * @related geom::Dual
 */
template <typename T, typename U, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> pow(U base, const geom::Dual<T,P> &xp) {
    T a_x = pow(base, xp.x);
    return geom::Dual<T>(
            a_x,
            xp.dx * a_x * log(base));
}


/** 
 * @brief Square root.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> sqrt(const geom::Dual<T,P> &d) {
    T sr = sqrt(d.x);
    return geom::Dual<T,P>(
            sr,
            d.dx / (2 * sr));
}


// min, max, ceil, floor

/**
 * @brief Absolute value.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> abs(const geom::Dual<T,P> &d) {
    bool neg = d.x < 0;
    T dx;
    
    if (d.x == 0) {
        dx = geom::Dual<T,P>::discontinuity_policy::resolve(-d.dx, d.dx);
    } else {
        dx = neg ? -d.dx: d.dx;
    }
    
    return geom::Dual<T,P>(neg ? -d.x : d.x, dx);
}


/**
 * @brief Ceiling function.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> ceil(const geom::Dual<T,P> &d) {
    T x = ceil(d.x);
    T dx;
    if (P == geom::DiscontinuityPolicy::Inf and x == d.x) {
        const T inf = std::numeric_limits<T>::infinity();
        // jump is in the same direction that the function is changing
        dx = std::copysign(inf, d.dx);
    } else if (P == geom::DiscontinuityPolicy::NaN and x == d.x) {
        dx = std::numeric_limits<T>::quiet_NaN();
    } else {
        dx = 0;
    }
    return geom::Dual<T,P>(x, dx);
}


/**
 * @brief Floor function.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> floor(const geom::Dual<T,P> &d) {
    T x = floor(d.x);
    T dx;
    if (P == geom::DiscontinuityPolicy::Inf and x == d.x) {
        const T inf = std::numeric_limits<T>::infinity();
        // jump is in the same direction that the function is changing
        dx = std::copysign(inf, d.dx);
    } else if (P == geom::DiscontinuityPolicy::NaN and x == d.x) {
        dx = std::numeric_limits<T>::quiet_NaN();
    } else {
        dx = 0;
    }
    return geom::Dual<T,P>(x, dx);
}


/**
 * @brief Minimum.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> min(const geom::Dual<T,P> &d1, const geom::Dual<T,P> &d2) {
    bool lt = d1.x <= d2.x;
    T x = lt ? d1.x : d2.x;
    T dx;
    if (d1.x == d2.x) {
        dx = geom::Dual<T,P>::discontinuity_policy::resolve(
            // faster increasing function wins min() on the left
            std::max(d1.dx, d2.dx),
            std::min(d1.dx, d2.dx));
    } else {
        dx = lt ? d1.dx : d2.dx;
    }
    return geom::Dual<T,P>(x,dx);
}


/**
 * @brief Maximum.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> max(const geom::Dual<T,P> &d1, const geom::Dual<T,P> &d2) {
    bool gt = d1.x >= d2.x;
    T x = gt ? d1.x : d2.x;
    T dx;
    if (d1.x == d2.x) {
        dx = geom::Dual<T,P>::discontinuity_policy::resolve(
            // slower increasing function wins max() on the left
            std::min(d1.dx, d2.dx),
            std::max(d1.dx, d2.dx));
    } else {
        dx = gt ? d1.dx : d2.dx;
    }
    return geom::Dual<T>(x,dx);
}


/**
 * @brief Fused multiply-add.
 * 
 * Compute `a * b + c`.
 * 
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> fma(
        geom::Dual<T,P> a,
        geom::Dual<T,P> b,
        geom::Dual<T,P> c)
{
    return geom::Dual<T,P>(
        geom::multiply_add(a.x, b.x, c.x),
        geom::sum_of_products(a.x, b.dx, a.dx, b.x) + c.dx);
}

// xxx: todo: implement proper decay semantics
template <typename T, typename U, geom::DiscontinuityPolicy P>
struct common_type<geom::Dual<T,P>,U> {
    typedef geom::Dual<typename common_type<T,U>::type, P> type;
};


template <typename T, geom::DiscontinuityPolicy P>
struct hash<geom::Dual<T,P>> {
    size_t operator()(const geom::Dual<T,P> &d) const {
        return geom::hash_combine(
            hash<T>{}(d.x),
            hash<T>{}(d.dx)
        );
    }
};

template <typename T, geom::DiscontinuityPolicy P>
inline bool isfinite(const geom::Dual<T,P> &d) {
    return std::isfinite(d.x) and std::isfinite(d.dx);
}

template <typename T, geom::DiscontinuityPolicy P>
inline bool isinf(const geom::Dual<T,P> &d) {
    return std::isinf(d.x) or std::isinf(d.dx);
}

template <typename T, geom::DiscontinuityPolicy P>
inline bool isnan(const geom::Dual<T,P> &d) {
    return std::isnan(d.x) or std::isnan(d.dx);
}

} // namespace std


#ifdef GEOMC_RANDOM_H_
#include <geomc/function/functiondetail/RandomDual.h>
#endif

// to avoid a circular include, but keep the functions organized correctly:
#include <geomc/function/Utils.h>

#endif	/* GEOMC_DUAL_H */

