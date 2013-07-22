/*
 *  Path.h
 * 
 * Represents a spline curve through N dimensions.
 * Knots have a position and a tangent.
 *
 *  Created on: Jul 14, 2011
 *      Author: tbabb
 */

#ifndef PATH_H_
#define PATH_H_

#include <vector>
#include <math.h>
#include "linalg/Vec.h"
#include "linalg/Ray.h"
#include "Expr.h"

//todo: A path could be a subclass of Raster with ptype=Ray!
//      ...except that Ray*float = Vec. Raster would assume Ray*float=Ray.
//      ...and a raster has fixed extent, while Path should allow dynamic resizing.
//      - Refactor raster to use decltype (I*O), i.e. allow distinct vertex point types and evaluation types.
//        > actually, what you want is <auto>, which is a c++11 only feature.
//        > Utils.h/interp_*() must use auto (or be return-type templated) too
//      - Refactor raster to allow dynamic storage; at least for 1d case.
//      Doing this will allow for modular interpolation schemes.
//      
//todo: option to parameterize by arc length
//      if A(s) is the arclength of the curve P(s),
//      then a length-parameterized curve is: P(A^-1(s))
//      if t(s) is the tangent vector, arclength is integral(sqrt(t(s)_x^2 + t(s)_y^2 + t(s)_z^2 + ... ), ds)
//      (in other words, the integral of the length of the tangent vector). 
//      There are suggestions that this is, in general, not integrable.
//      For a linear-interpolated spline curve, this reduces to the integral of the square root of a
//      quadratic function, in other words a hyperbola. The hyperbola happens to be positive and real-valued.
//      according to wolframpalpha, the solution seems to be of the form
//      ||S_tangent|| * (0.5*s + a) + b*sinh^-1(cx +ds). The formula for a, b, c, and d are not given.
//      In the general case, the computed indefinite integral is much uglier.
//todo: option to split at a particular point (in both params)
//todo: option to split tangents at a particular knot
//todo: algorithm to rasterize by subdividing recursively into line segments.
//todo: capability to change ease fn? (linear ease == quadratic/s^2 overall)
//todo: evaluate tangent? (return an Expr?)
//      (not too hard; derivative of P(s))
//todo: evaluate acceleration? (as Expr?)
//todo: move this out of Function and into Shape.
//todo: allow sensible 1D usage 
//todo: with a width spline, trace()
//todo: bound() (fits in the convex hull of control points)
//      could solve for X and Y extremes?

namespace geom {

template <typename T, index_t N>
class Path {
public:
    Path(){}

    Vec<T,N> eval(T t){
        size_t n = t;
        size_t sz = knots.size();
        
        //degenerate cases
        if (sz < 1){
            return Vec<T,N>();
        } else if (t < 0){ 
            return knots[0].atMultiple(t);
        } else if (n >= sz - 1){
            T t_0 = t - sz + 1;
            return knots[sz - 1].atMultiple(t_0);
        }
        
        T t_0 = t - n;
        Ray<T,N> r0 = knots[n];
        Ray<T,N> r1 = knots[n+1];
        Vec<T,N> v0 = r0.atMultiple(t_0);
        Vec<T,N> v1 = r1.atMultiple(t_0 - 1);
        return v0.mix(v1, ease(t_0));
    }
    
    void appendKnot(Vec<T,N> p, Vec<T,N> v){
        knots.push_back(Ray<T,N>(p,v));
    }
    
    std::vector< Ray<T,N> > knots;
    
protected:
    static inline T ease(T t){
        return t*t*t*(t*(t*6 - 15) + 10);
    }
};

} /* end namespace geom */

#endif /* PATH_H_ */
