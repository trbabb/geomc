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
 * @brief A wrapper shape which orients another arbitrary shape with an AffineTransform.
 * 
 * Oriented shapes can be constructed simply by applying an AffineTransform to an
 * ordinary shape:
 *
 *     AffineTransform<double,3> xf = rotation(...);
 *     Oriented<Cylinder<double,3>> oriented = xf * Cylinder<double,3>();
 *
 * Transforming an Oriented results in another Oriented of the same type:
 *
 *     Oriented<Cylinder<double,3>> ocyl_a = xf * Cylinder<double,3>();
 *     Oriented<Cylinder<double,3>> ocyl_b = xf * ocyl_a;
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
    
    Rect<T,N> bounds() const {
        // the base class implementation would work; but we specialize
        // here for a little extra performance (avoid stupidly re-computing `d`).
        Rect<T,N> r;
        for (index_t axis = 0; axis < N; ++axis) {
            // we want to test along a world cardinal axis in shape-space coordinates.
            // but this direction is a normal, which is transformed
            // with the inverse transpose. so we take a row of xf.mat instead of 
            // a col from xf.inv:
            Vec<T,N> basis(xf.mat.row(axis));
            // find the body-space extreme and transform it to world space:
            Vec<T,N> p_hi = xf * shape.convex_support( basis);
            Vec<T,N> p_lo = xf * shape.convex_support(-basis);
            // update the extremes along the test axis
            r.hi[axis]  = p_hi[axis];
            r.lo[axis]  = p_lo[axis];
        }
        return r;
    }
    
    /// Return an Oriented Rect containing the shape.
    inline Oriented<Rect<T,N>> oriented_bounds() const {
        return xf * shape.bounds();
    }
    
    
    /// Transform this `Oriented` by `xf`
    Oriented<Shape>& operator*=(const AffineTransform<T,N>& xf) {
        this->xf *= xf;
        return *this;
    }
    
    
    /// Transform this `Oriented` by the inverse of `xf`.
    Oriented<Shape>& operator/=(const AffineTransform<T,N>& xf) {
        this->xf /= xf;
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