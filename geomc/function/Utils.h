/*
 * Utils.h
 *
 *  Created on: Feb 22, 2009
 *      Author: Tim Babb
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <cmath>
#include <algorithm>

#include <geomc/function/functiondetail/UtilDetail.h>

namespace geom {
    
    /**
     * @ingroup function
     * @{
     */
    
    /**
     * A high-precision method for computing (a * b) - (c * d). In cases
     * where the two products are large and close in value, the result can be
     * wildly wrong due to catastrophic cancellation when using the "naive" method.
     * This function will return a result within 1.5 ULP of the exact answer.
     *
     * A high-precision method is available for float, double, and long double.
     * Other types will fall back to the "naive" method.
     */
    template <typename T> 
    inline T diff_of_products(T a, T b, T c, T d) {
        return detail::_ImplFMADot<T>::diff_of_products(a,b,c,d);
    }
    
    /**
     * A high-precision method for computing (a * b) + (c * d). In cases
     * where the two products are large and similar in magnitude but opposite in sign, 
     * the result can be wildly wrong due to catastrophic cancellation when using 
     * the "naive" method. This function will return a result within 1.5 ULP of the 
     * exact answer.
     *
     * A high-precision method is available for float, double, and long double.
     * Other types will fall back to the "naive" method.
     */
    template <typename T>
    inline T sum_of_products(T a, T b, T c, T d) {
        return detail::_ImplFMADot<T>::sum_of_products(a,b,c,d);
    }

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
    inline T mix(S t, T a, T b) {
        return sum_of_products(1 - t, a, t, b);
    }
    
    /**
     * Linearly interpolate across `dim` input dimensions.
     * @param s An array of `dim` values; a point on `(0, 1)`<sup>`dim`</sup> at which to interpolate.
     * @param p An array of the 2<sup>`dim`</sup> values adjacent to `s`.
     * @param dim Number of input dimensions.
     * @return Interpolated value.
     */
    template <typename S, typename T>
    T interp_linear(const S s[], const T p[], int dim) {
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
     * @param s Interpolation parameter between 0 and 1, representing the point between p[1] and p[2] at which to interpolate.
     * @param p An array of 4 values surrounding `s`.
     * @return Interpolated value.
     */
    template <typename S, typename T>
    T interp_cubic(S s, const T p[4]) {
        return p[1] + 0.5 * s*(p[2] - p[0] + s*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + s*(3.0*(p[1] - p[2]) + p[3] - p[0])));
    }
    
    /**
     * Cubic interpolation across `dim` input dimensions.
     * @param s An array of `dim` values; a point on `(0, 1)`<sup>`dim`</sup> representing the position in the center grid cell at which to interpolate.
     * @param p An array of the 4<sup>dim</sup> values adjacent to `s`.
     * @param dim Number of input dimensions.
     * @return Interpolated value.
     */
    template <typename S, typename T>
    T interp_cubic(const S s[], const T p[], int dim) {
        if (dim <= 1){
            return interp_cubic(*s, p);
        } else {
            T buf[4];
            int next_dim = dim - 1;
            int blocksize = 1 << (2*next_dim); //TODO: is this right?
            buf[0] = interp_cubic(s, p,               next_dim);
            buf[1] = interp_cubic(s, p +   blocksize, next_dim);
            buf[2] = interp_cubic(s, p + 2*blocksize, next_dim);
            buf[3] = interp_cubic(s, p + 3*blocksize, next_dim);
            return interp_cubic(s[next_dim], buf);
        }
    }
    
    /**
     * Solve the quadratic equation `ax`<sup>2</sup>` + bx + c  = 0`.
     * Formulated to preserve precision; the "high school" formula is unstable.
     * 
     * @param [out] results Roots of the specified quadratic equation.
     * @param a Coefficient of `x`<sup>2</sup>.
     * @param b Coefficient of `x`.
     * @param c Constant.
     * @return `true` if real roots exist; `false` otherwise (in which case 
     * the contents of `results` will be unaltered).
     */
    template <typename T> bool quadratic_solve(T results[2], T a, T b, T c){
        T descr = diff_of_products(b, b, 4 * a, c);
        if (descr < 0) return false; //no real solution
        T q = -0.5 * (b + copysign(sqrt(descr), b));
        results[0] = q / a;
        results[1] = c / q;
        return true;
    }
    
    /// @} // ingroup function
}

#endif /* UTILS_H_ */
