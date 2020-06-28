/* 
 * File:   Oriented.h
 * Author: tbabb
 */

#ifndef ORIENTED_H
#define ORIENTED_H

#include <geomc/linalg/AffineTransform.h>
#include <geomc/shape/Bounded.h>

namespace geom {
    
// todo: use std::forward<X>() on constructors, to allow move-construction
// todo: use ::N and ::elem_t with Extrusion to make things easier and
//       more consistent
// todo: test that the above scheme is actually legal c++11
// todo: specializations for OrientedRect
//       > incl. transform fns
// todo: implement sdf()
// todo: implement nearest(), which admits sdf() for compositional shapes


template <typename Shape>
class Oriented: public virtual Convex<typename Shape::elem_t, Shape::N> {
public:
    typedef typename Shape::elem_t T;
    static constexpr size_t N = Shape::N;
    
    AffineTransform<T,N> xf;
    Shape shape;
    
    
    Oriented(const Shape& s):shape(s) {}
    
    Oriented(const Shape& s, const AffineTransform<T,N>& xf):
        xf(xf),
        shape(s) {}
        
    
    bool contains(Vec<T,N> p) const {
        p = p / xf;
        return shape.contains(p);
    }
    
    
    Vec<T,N> convex_support(Vec<T,N> d) const {
        d = xf.applyInverseNormal(d);
        Vec<T,N> p = shape.convex_support(d);
        return xf * p;
    }
    
    
    Oriented<Shape>& operator*=(const AffineTransform<T,N>& xf1) {
        xf *= xf1;
        return *this;
    }
    
    
    Oriented<Shape>& operator/=(const AffineTransform<T,N>& xf1) {
        xf /= xf1;
        return *this;
    }
    
}; // class Oriented


template <typename Shape>
Oriented<Shape> operator*(
        const AffineTransform<typename Shape::elem_t, Shape::N>& xf,
        const Shape& s) {
    return Oriented<Shape>(s);
}


template <typename Shape>
Oriented<Shape> operator*(
        const AffineTransform<typename Shape::elem_t, Shape::N>& xf,
        const Oriented<Shape>& s) {
    return Oriented<Shape>(s, xf * s.xf);
}


} // namespace geom

#endif // ORIENTED_H