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

namespace geom {

    template <typename T>
    inline T clamp(const T &v, const T &lo, const T &hi){
        return std::max(lo, std::min(hi, v));
    }
    
    template <typename S, typename T>
    inline T mix(S t, T a, T b){
        return (1-t)*a + t*b;
    }
    
    template <typename S, typename T>
    T interp_linear(const S s[], const T p[], int dim){
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
    
    template <typename S, typename T>
    T interp_cubic(S s, const T p[4]){
        return p[1] + 0.5 * s*(p[2] - p[0] + s*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + s*(3.0*(p[1] - p[2]) + p[3] - p[0])));
    }
    
    template <typename S, typename T>
    T interp_cubic(const S s[], const T p[], int dim){
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
    
    // formulated to preserve precision; "high school" formula is unstable.
    template <typename T> bool quadratic_solve(T results[2], T a, T b, T c){
        T descr = b*b - 4*a*c;
        if (descr < 0) return false; //no real solution
        T q = -0.5 * (b + copysign(sqrt(descr), b));
        results[0] = q/a;
        results[1] = c/q;
        return true;
    }
}

#endif /* UTILS_H_ */
