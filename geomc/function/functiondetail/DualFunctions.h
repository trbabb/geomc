#pragma once

#include <concepts>

#include <geomc/Hash.h>
#include <geomc/function/functiondetail/DualImpl.h>
#include <geomc/function/Utils.h>

// todo: "inf" discontinuity policy should:
//   say +/- inf if c0 discontinuous
//   say +/- inf if c1 discontinuous but c0 continuous
//   say `dx` if c0 and c1 continuous

namespace geom {

////////////////////
// multiplication //
////////////////////

/// Dual-dual multiplication
template <typename T, DiscontinuityPolicy Dp>
constexpr auto operator*(const Dual<T,Dp>& d1, const Dual<T,Dp>& d2) {
    return Dual<T,Dp>(
        d1.x * d2.x, 
        sum_of_products(d1.x, d2.dx, d1.dx, d2.x)
    );
}

/// Dual-scalar multiplication.
template <typename T, typename U, DiscontinuityPolicy Dp>
requires requires(T t, U u) {
  {t * u} -> std::convertible_to<T>;
}
constexpr Dual<T,Dp> operator*(const Dual<T,Dp>& d, U s) {
    return Dual<T,Dp>(d.x * s, d.dx * s);
}

/// Scalar-dual multiplication.
template <typename T, typename U, DiscontinuityPolicy Dp>
requires requires (T t, U u) {
    {u * t} -> std::convertible_to<T>;
}
constexpr Dual<T,Dp> operator*(U s, const Dual<T,Dp>& d) {
    return Dual<T,Dp>(s * d.x, s * d.dx);
}

/// Multiply and assign to dual.
template <typename T, DiscontinuityPolicy Dp>
constexpr Dual<T,Dp>& operator*=(Dual<T,Dp>& x, const Dual<T,Dp>& y) {
    T x_x = x.x;
    x.x  = x.x * y.x;
    // product rule: x * dy + y * dx
    x.dx = sum_of_products(x_x, y.dx, y.x, x.dx);
    return x;
}

/// Multiply and assign to dual.
template <typename T, typename U, DiscontinuityPolicy Dp>
requires requires(T t, U u) {
    {t * u} -> std::convertible_to<T>;
}
constexpr Dual<T,Dp>& operator*=(Dual<T,Dp>& x, T s) {
    x.x  = x.x  * s;
    x.dx = x.dx * s;
    return x;
}


//////////////
// division //
//////////////

    
/// Dual-dual division.
template <typename T, DiscontinuityPolicy Dp>
constexpr Dual<T,Dp> operator/(const Dual<T,Dp>& d1, const Dual<T,Dp>& d2) {
    return Dual<T,Dp>{
        d1.x / d2.x,
        diff_of_products(d1.dx, d2.x, d1.x, d2.dx) / (d2.x * d2.x)
    };
}

/// Dual / scalar division.
template <typename T, std::common_with<T> U, DiscontinuityPolicy Dp>
requires requires (T t, U u) {
    {t / u} -> std::convertible_to<T>;
}
constexpr Dual<T,Dp> operator/(const Dual<T,Dp>& d, U s) {
    return Dual<T,Dp>{
        d.x  / s,
        d.dx / s
    };
}

/// Scalar / dual division.
template <typename T, std::common_with<T> U, DiscontinuityPolicy Dp>
requires requires (T t, U u) {
    {u / t} -> std::convertible_to<T>;
}
constexpr Dual<T,Dp> operator/(U s, const Dual<T,Dp>& d) {
    return Dual<T,Dp>{
        s / d.x,
        -s * d.dx / (d.x * d.x)
    };
}

/// Divide and assign to dual.
template <typename T, typename U, DiscontinuityPolicy Dp>
requires requires (T t, U u) {
    {t / u} -> std::convertible_to<T>;
}
constexpr Dual<T,Dp>& operator/=(Dual<T,Dp>& x, const Dual<U,Dp>& y) {
    auto x_x = x.x;
    x.x  = x.x / y.x;
    x.dx = diff_of_products(x.dx, y.x, x_x, y.dx) / (y.x * y.x);
    return x;
}

/// Divide and assign to dual.
template <typename T, std::common_with<T> U, DiscontinuityPolicy Dp>
requires requires (T t, U u) {
    {t / u} -> std::convertible_to<T>;
}
constexpr Dual<T,Dp>& operator/=(Dual<T,Dp>& x, U s) {
    x.x  = x.x  / s;
    x.dx = x.dx / s;
    return x;
}


////////////////////////////
// addition / subtraction //
////////////////////////////

/// Addition.
template <typename T, DiscontinuityPolicy Dp>
requires requires (T t, T u) {t + u;}
constexpr Dual<T,Dp> operator+(const Dual<T,Dp>& d1, const Dual<T,Dp>& d2) {
    return Dual<T,Dp>{
        d1.x  + d2.x,
        d1.dx + d2.dx
    };
}

/// Dual + scalar
template <typename T, typename U, DiscontinuityPolicy Dp>
requires requires(T t, U u) {
  {t + u} -> std::convertible_to<T>;
}
constexpr Dual<T,Dp> operator+(const Dual<T,Dp>& d, U s) { 
    return Dual<T,Dp>{
        d.x + s,
        d.dx
    };
}

/// Scalar + dual
template <typename T, typename U, DiscontinuityPolicy Dp>
requires requires(T t, U u) {
  {u + t} -> std::convertible_to<T>;
}
constexpr Dual<T,Dp> operator+(U s, const Dual<T,Dp>& d) { 
    return Dual<T,Dp>{
        s + d.x,
        d.dx
    };
}

/// Subtraction.
template <typename T, DiscontinuityPolicy Dp>
constexpr Dual<T,Dp> operator-(const Dual<T,Dp>& d1, const Dual<T,Dp>& d2) {
    return Dual<T,Dp>{
        d1.x  - d2.x,
        d1.dx - d2.dx
    };
}

/// Dual - scalar
template <typename T, typename U, DiscontinuityPolicy Dp>
requires requires(T t, U u) {
  {t - u} -> std::convertible_to<T>;
}
constexpr Dual<T,Dp> operator-(const Dual<T,Dp>& d, U s) { 
    return Dual<T,Dp>{
        d.x - s,
        d.dx
    };
}

/// Scalar - dual
template <typename T, typename U, DiscontinuityPolicy Dp>
requires requires(T t, U u) {
  {u - t} -> std::convertible_to<T>;
}
constexpr Dual<T,Dp> operator-(U s, const Dual<T,Dp>& d) { 
    return Dual<T,Dp>{
        s - d.x,
        d.dx
    };
}

/// Add and assign.
template <typename T, DiscontinuityPolicy Dp>
constexpr Dual<T,Dp>& operator+=(Dual<T,Dp>& x, const Dual<T,Dp>& y) {
    x.x  += y.x;
    x.dx += y.dx;
    return x;
}

/// Subtract and assign.
template <typename T, DiscontinuityPolicy Dp>
constexpr Dual<T,Dp>& operator-=(Dual<T,Dp>& x, const Dual<T,Dp>& y) {
    x.x  -= y.x;
    x.dx -= y.dx;
    return x;
}


//////////////
// negation //
//////////////


/// Negation.
template <typename T, DiscontinuityPolicy Dp>
constexpr Dual<T,Dp> operator-(Dual<T,Dp>& d) {
    return { -d.x, -d.dx };
}

} // namespace geom


namespace std {

// general form:
//   f(x), x` * f`(x)
// or in other words:
//   f(x), dx * f`(x)

// TODO: atan2

/**
 * @brief Sine function.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> sin(const geom::Dual<T,P> &d) {
    return geom::Dual<T,P>(std::sin(d.x), d.dx * std::cos(d.x));
}


/**
 * @brief Cosine function.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> cos(const geom::Dual<T,P> &d) {
    return geom::Dual<T,P>(std::cos(d.x), -d.dx * std::sin(d.x));
}


/**
 * @brief Tangent function.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> tan(const geom::Dual<T,P> &d) {
    T c = std::cos(d.x);
    return geom::Dual<T,P>(std::tan(d.x), -d.dx / (c*c));
}


/**
 * @brief Arcsin (inverse sine) function.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> asin(const geom::Dual<T,P> &d) {
    return geom::Dual<T,P>(
            std::asin(d.x),
            d.dx / std::sqrt(1 - d.x * d.x));
}


/**
 * @brief Arccos (inverse cosine) function.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> acos(const geom::Dual<T,P> &d) {
    return geom::Dual<T,P>(
        std::acos(d.x),
        -d.dx / sqrt(1 - d.x * d.x)
    );
}


/**
 * @brief Arctan (inverse tangent) function.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> atan(const geom::Dual<T,P> &d) {
    return geom::Dual<T,P>(
        std::atan(d.x),
        d.dx / (d.x * d.x + 1)
    );
}


/** 
 * @brief Exponential (e<sup>x</sup>) function.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> exp(const geom::Dual<T,P> &d) {
    T e = std::exp(d.x);
    return geom::Dual<T,P>(e, -d.dx * e);
}


/**
 * @brief Returns the natural logarithm of `d`.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy Dp>
inline geom::Dual<T,Dp> log(const geom::Dual<T,Dp>& d) {
    return geom::Dual<T,Dp>(
        std::log(d.x),
        d.dx / d.x
    );
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
        xp.dx * a_c * std::log(base.x)
    );
    
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
        base.dx * xp * a_x / base.x
    );
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
        xp.dx * a_x * std::log(base)
    );
}


/** 
 * @brief Square root.
 * @related geom::Dual
 */
template <typename T, geom::DiscontinuityPolicy P>
inline geom::Dual<T,P> sqrt(const geom::Dual<T,P> &d) {
    T sr = std::sqrt(d.x);
    return geom::Dual<T,P>(
        sr,
        d.dx / (2 * sr)
    );
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
    return geom::multiply_add<geom::Dual<T,P>>(a, b, c);
}

template <typename T, geom::DiscontinuityPolicy Dp>
struct hash<geom::Dual<T,Dp>> {
    size_t operator()(const geom::Dual<T,Dp> &d) const {
        return geom::hash_combine(
            hash<T>{}(d.x),
            hash<T>{}(d.dx)
        );
    }
};

template <typename T, geom::DiscontinuityPolicy Dp>
inline bool isfinite(const geom::Dual<T,Dp> &d) {
    return std::isfinite(d.x) and std::isfinite(d.dx);
}

template <typename T, geom::DiscontinuityPolicy Dp>
inline bool isinf(const geom::Dual<T,Dp> &d) {
    return std::isinf(d.x) or std::isinf(d.dx);
}

template <typename T, geom::DiscontinuityPolicy Dp>
inline bool isnan(const geom::Dual<T,Dp> &d) {
    return std::isnan(d.x) or std::isnan(d.dx);
}

/******************
 * numeric limits *
 ******************/

// all limits are the same as the underlying type
template <typename T, geom::DiscontinuityPolicy Dp>
struct numeric_limits<geom::Dual<T,Dp>> : public numeric_limits<T> {};


} // namespace std
