#pragma once

/*
 * Utils.h
 *
 *  Created on: Feb 22, 2009
 *      Author: Tim Babb
 */

#include <cmath>
#include <algorithm>

#include <geomc/function/functiondetail/UtilDetail.h>

namespace geom {

/**
 * @ingroup function
 * @{
 */


template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, T> positive_mod(T a, T b) {
    return (a % b + b) % b;
}

template <typename T>
inline std::enable_if_t<not std::is_integral_v<T>, T> positive_mod(T a, T b)
{
    return a - b * std::floor(a / b);
}

/**
 * A high-precision method for computing (a * b) - (c * d). In cases
 * where the two products are large and close in value, the result can be
 * wildly wrong due to catastrophic cancellation when using the "naive" method.
 * This function will return a result within 1.5 ULP of the exact answer.
 *
 * A high-precision method is available for float, double, long double, and Dual numbers
 * composed of those types. Other types will fall back to the "naive" method.
 */
template <typename T> 
inline T diff_of_products(T a, T b, T c, T d) {
    return detail::_ImplFMADot<T>::diff_of_products(a, b, c, d);
}

/**
 * A high-precision method for computing (a * b) + (c * d). In cases
 * where the two products are large and similar in magnitude but opposite in sign, 
 * the result can be wildly wrong due to catastrophic cancellation when using 
 * the "naive" method. This function will return a result within 1.5 ULP of the 
 * exact answer.
 *
 * A high-precision method is available for float, double, long double, and Dual numbers
 * composed of those types. Other types will fall back to the "naive" method.
 */
template <typename T>
inline T sum_of_products(T a, T b, T c, T d) {
    return detail::_ImplFMADot<T>::sum_of_products(a, b, c, d);
}

/**
 * @brief Fused multiply-add.
 * 
 * Compute `a * b + c` using a fused multiply add operation (if available), which may be both more
 * performant and more precise than multiple instructions.
 * 
 * Overloaded such that the calculation falls back to an ordinary composite multiply-add if a fused
 * instruction is not available for `T`. (Compare to `std::fma()` which will promote to floating point,
 * possibly losing precision).
 */
template <typename T>
inline T multiply_add(T a, T b, T c) {
    return detail::_ImplFMADot<T>::multiply_add(a, b, c);
}

// todo: it would be more consistent with the shape sub-library
//  to call the following "clip":

/**
 * Clamp `v` between `lo` and `hi`.
 * @param v Value to clamp.
 * @param lo Minimum value of `v`.
 * @param hi Maximum value of `v`.
 */
template <typename T>
inline T clamp(const T &v, const T &lo, const T &hi) {
    return std::max(lo, std::min(hi, v));
}

/**
 * Linearly interpolate between `a` and `b`.
 * @param t Interpolation parameter, between 0 and 1.
 * @param a Value at `t = 0`.
 * @param b Value at `t = 1`.
 * @return Interpolated value.
 */
template <typename S, typename T>
inline typename std::common_type<S,T>::type mix(S t, T a, T b) {
    typedef typename std::common_type<S,T>::type U;
    // (1 - t) * a + t * b
    return sum_of_products<U>(1 - t, a, t, b);
}

/**
 * Linearly interpolate across `dim` input dimensions.
 * @param s An array of `dim` values; a point on `(0, 1)`<sup>`dim`</sup> at which to
 * interpolate.
 * @param p An array of the 2<sup>`dim`</sup> values adjacent to `s`.
 * @param dim Number of input dimensions.
 * @return Interpolated value.
 */
template <typename S, typename T>
inline T interp_linear(const S s[], const T p[], int dim) {
    if (dim <= 1){
        return mix<S,T>(*s, p[0], p[1]);
    } else {
        int next_dim = dim - 1;
        int blocksize = 1 << next_dim;
        T p0 = interp_linear(s, p,             next_dim);
        T p1 = interp_linear(s, p + blocksize, next_dim);
        return mix<S,T>(s[next_dim], p0, p1);
    }
}

/**
 * Cubic interpolation in 1 dimension.
 * @param s Interpolation parameter between 0 and 1, representing the point
 * between p[1] and p[2] at which to interpolate.
 * @param p An array of 4 values surrounding `s`.
 * @return Interpolated value.
 */
template <typename S, typename T>
inline T interp_cubic(S s, const T p[4]) {
    return p[1] + 0.5 * 
        s*(p[2] - p[0] + 
            s*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + 
                s*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}

/**
 * Cubic interpolation across `dim` input dimensions.
 * @param s An array of `dim` values; a point on `(0, 1)`<sup>`dim`</sup> representing 
 * the position in the center grid cell at which to interpolate.
 * @param p An array of the 4<sup>dim</sup> values adjacent to `s`.
 * @param dim Number of input dimensions.
 * @return Interpolated value.
 */
template <typename S, typename T>
inline T interp_cubic(const S s[], const T p[], int dim) {
    if (dim <= 1){
        return interp_cubic(*s, p);
    } else {
        T buf[4];
        int next_dim = dim - 1;
        int blocksize = 1 << (2 * next_dim); //TODO: is this right?
        buf[0] = interp_cubic(s, p,                 next_dim);
        buf[1] = interp_cubic(s, p +     blocksize, next_dim);
        buf[2] = interp_cubic(s, p + 2 * blocksize, next_dim);
        buf[3] = interp_cubic(s, p + 3 * blocksize, next_dim);
        return interp_cubic(s[next_dim], buf);
    }
}

/**
 * Solve the quadratic equation `ax`<sup>2</sup>` + bx + c  = 0`.
 * Formulated to preserve precision and avoid any divisions by zero.
 * 
 * In the case where there is exactly one real solution, the two roots
 * will be set equal to each other.
 * 
 * @param [out] results Roots of the specified quadratic equation.
 * @param a Coefficient of `x`<sup>2</sup>.
 * @param b Coefficient of `x`.
 * @param c Constant.
 * @return `true` if real roots exist and the coefficients are not all zero;
 * `false` otherwise (in which case the contents of `results` will be unaltered).
 */
template <typename T>
inline bool quadratic_solve(T results[2], T a, T b, T c){
    T descr = diff_of_products(b, b, 4 * a, c); // b^2 - 4ac
    // no real solution, or an equation of the form "c = 0"
    if (descr < 0 or (a == 0 and b == 0)) return false;
    // addition cannot cancel, since b and descr are chosen to have the same sign
    T q = -0.5 * (b + std::copysign(std::sqrt(descr), b));
    T x0 = q / a; // note: client could pass a = 0
    T x1 = c / q; // q is zero iff b and c are zero
    results[0] = (a != 0) ? x0 : x1;
    results[1] = (q != 0) ? x1 : x0;
    return true;
}

/**
 * Compute the shortest angle from `a0` to `a1`, in radians.
 * 
 * Returns a value in the range `[-π, π]`.
 */
template <typename T>
inline T angle_to(T radians_0, T radians_1) {
    constexpr T turn = 2 * M_PI;
    T a0 = positive_mod<T>(radians_0, turn);
    T a1 = positive_mod<T>(radians_1, turn);
    // one of three images of a0 is closest to a1:
    // ···----o------|----o-->*----|---o---···
    //   a0 - τ      0   a0   a1   τ   a0 + τ  
    T da = a1 -  a0;
    T db = a1 - (a0 + turn);
    T dc = a1 - (a0 - turn);
    // find the difference of smallest magnitude:
    T err = std::abs(da)  < std::abs(db) ? da  : db;
    err   = std::abs(err) < std::abs(dc) ? err : dc;
    return err;
}


/**
 * @brief Compute `ceil(a/b)`.
 * 
 * @param a Numerator.
 * @param b Denominator.
 */
template <typename T>
constexpr auto ceil_div(T a, T b) {
    if constexpr (not std::is_integral<T>()) {
        return std::ceil(a / b);
    } else {
        return (a + b - 1) / b;
    }
}

/**
 * @brief Compute `floor(a/b)`.
 * 
 * @param a Numerator.
 * @param b Denominator.
 */
template <typename T>
constexpr auto floor_div(T a, T b) {
    if constexpr (not std::is_integral<T>()) {
        return std::floor(a / b);
    } else {
        return (a + b) / b - 1;
    }
}

template <typename T>
constexpr bool is_power_of_two(T x) requires std::integral<T> {
    return x > 0 and (x & (x - 1)) == 0;
}

/// @} // ingroup function

} // namespace geom
