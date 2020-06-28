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
// todo: specializations for OrientedRect
//       > incl. transform fns
// todo: implement sdf()
// todo: implement nearest(), which admits sdf() for compositional shapes

/**
 * @ingroup shape
 * @brief A wrapper shape which orients another shape with an AffineTransform.
 * 
 * Oriented shapes can be constructed simply by applying an AffineTransform to an
 * ordinary shape:
 *
 *     AffineTransform<double,3> xf = rotation(...);
 *     Oriented<Cylinder<double,3>> oriented = xf * Cylinder<double,3>();
 *
 * Transforming an Oriented results in another Oriented of the same type:
 *
 *     Oriented<Cylinder<double,3>> ocyl1 = xf * Cylinder<double,3>();
 *     Oriented<Cylinder<double,3>> ocyl2 = xf * ocyl1;
 * 
 */
template <typename Shape>
class Oriented: public virtual Convex<typename Shape::elem_t, Shape::N> {
public:
    /// The coordinate type of this Shape
    typedef typename Shape::elem_t T;
    /// The dimension of this Shape
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


/**
 * @brief Transform the shape `s` by wrapping it in an `Oriented` class.
 * @related AffineTransform
 * @related Oriented
 */
template <typename Shape>
Oriented<Shape> operator*(
        const AffineTransform<typename Shape::elem_t, Shape::N>& xf,
        const Shape& s) {
    return Oriented<Shape>(s, xf);
}


/**
 * @brief Transform the oriented shape `s` by `xf`.
 * @related AffineTransform
 * @related Oriented
 */
template <typename Shape>
Oriented<Shape> operator*(
        const AffineTransform<typename Shape::elem_t, Shape::N>& xf,
        const Oriented<Shape>& s) {
    return Oriented<Shape>(s, xf * s.xf);
}


} // namespace geom

#endif // ORIENTED_H