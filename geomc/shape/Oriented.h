#pragma once

#include <geomc/linalg/AffineTransform.h>
#include <geomc/shape/Rect.h>
#include <geomc/shape/shapedetail/SeparatingAxis.h>


namespace geom {

// todo: make Rotor, and with it make SimilarTransform (Rotor + Translate), 
//    and have Oriented store a SimilarTransform instead of an AffineXf
//    > possibly internally delegate to Transformed<...> for efficiency?
//    > allow sdf() and project()
// todo: rename current Oriented impl to Transformed
// todo: implement sdf() on Rotor-based Oriented
//   - (can't work for matrix-based due to shears and scales)
//   - (until then, user can do it themselves if they know there is no scale)
// todo: implement project() on Rotor-based Oriented
//   - (can't work for matrix-based due to shears and scales)
// todo: use std::forward<X>() on constructors, to allow move-construction

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
class Oriented:
    public Convex          <typename Shape::elem_t, Shape::N, Oriented<Shape>>,
    public RayIntersectable<typename Shape::elem_t, Shape::N, Oriented<Shape>>
{
public:
    /// The coordinate type of this Shape
    typedef typename Shape::elem_t T;
    /// The dimension of this Shape
    static constexpr size_t N = Shape::N;
    
    /// Un-transformed shape.
    Shape shape;
    /// Transformation orienting `shape`.
    AffineTransform<T,N> xf;
    
    /// Wrap a default-constructed Shape with the identity transformation.
    Oriented() {}
    
    /// Wrap `s` in an identity transformation.
    explicit Oriented(const Shape& s):shape(s) {}
    
    /// Construct an oriented shape from `s`, which is positioned and oriented by `xf`.
    Oriented(const Shape& s, const AffineTransform<T,N>& xf):
        shape(s),
        xf(xf) {}
    
    bool operator==(const Oriented& other) const {
        return shape == other.shape && xf == other.xf;
    }
    
    /// Shape-point overlap test.
    bool contains(Vec<T,N> p) const {
        p = p / xf;
        return shape.contains(p);
    }
    
    Vec<T,N> convex_support(Vec<T,N> d) const {
        d = xf.apply_inverse_normal(d);
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
    
    /// Ray-shape intersection.
    Rect<T,1> intersect(const Ray<T,N>& r) const {
        return shape.intersect(r / xf);
    }
    
    /// Return an Oriented Rect containing the shape.
    inline Oriented<Rect<T,N>> oriented_bounds() const {
        return xf * shape.bounds();
    }
    
    
}; // class Oriented


/**
 * @ingroup shape
 * @brief Partial specialization of Oriented for Rects.
 */
template <typename T, index_t N>
class Oriented<Rect<T,N>>:
    public Convex<T,N,Oriented<Rect<T,N>>>,
    public RayIntersectable<T,N,Oriented<Rect<T,N>>>
{    
public:
    
    /// Un-transformed extents.
    Rect<T,N> shape;
    /// Transformation orienting `shape`.
    AffineTransform<T,N> xf; // moves object points to world positions
    
    /// Construct an empty axis-aligned box.
    Oriented() {}
    
    /// Construct an axis-aligned box from the given Rect.
    explicit Oriented(const Rect<T,N>& box):shape(box) {}
    
    /// Construct an oriented box from the given Rect and object-to-world transformation.
    Oriented(const Rect<T,N>& box, const AffineTransform<T,N>& xf):
            shape(box),
            xf(xf) {}
    
    bool operator==(const Oriented& other) const {
        return shape == other.shape && xf == other.xf;
    }
    
    /// Obtain an axis-aligned bounding box for this shape.
    Rect<T,N> bounds() const {
        // faster than using convex_support, somewhat surprisingly.
        Vec<T,N> pts[2] = {shape.lo, shape.hi};
        Vec<T,N> lo = xf * shape.lo;
        Vec<T,N> hi = lo;
        
        // for each corner
        for (index_t i = 1; i < (1 << N); i++) {
            Vec<T,N> p;
            // for each coordinate
            for (index_t j = 0; j < N; j++) {
                p[j] = pts[(i & (1 << j)) != 0][j];
            }
            p  = xf * p;
            lo = std::min(lo, p);
            hi = std::max(hi, p);
        }
        
        return Rect<T,N>(lo, hi);
    }

    /// Return `this`.
    inline Oriented< Rect<T,N> > oriented_bounds() const {
        return *this;
    }
    
    Vec<T,N> convex_support(Vec<T,N> d) const {
        Vec<T,N> d_body = xf.apply_inverse_normal(d);
        Vec<T,N> o;
        for (index_t i = 0; i < N; i++) {
            o[i] = ((d_body[i] < 0) ? shape.lo : shape.hi)[i];
        }
        return xf * o;
    }
    
    /**
     * @brief Returns true if and only if `p` is inside or on the surface of 
     * this shape.
     */
    bool contains(Vec<T,N> p) const {
        p = p / xf; // world coord to object coord
        return shape.contains(p);
    }
    
    /**
     * Test whether this Oriented<Rect> overlaps an axis-aligned Rect.
     * 
     * @param b1 OritentedRect to test against.
     * @return `true` if and only if `this` overlaps with `b1`; `false` otherwise.
     */
    bool intersects(const Rect<T,N>& r) const {
        return detail::RectIntersector<T,N>::intersect(*this, r);
    }
    
    /**
     * Test whether this Oriented<Rect> overlaps another.
     * 
     * @param b1 OritentedRect to test against.
     * @return `true` if and only if `this` overlaps with `b1`; `false` otherwise.
     */
    bool intersects(const Oriented< Rect<T,N> >& b1) {
        // we will make ourselves axis-aligned; this is the same as applying
        // the inverse of our xf. We apply the inverse of xf to b1 too,
        // to preserve the relationhip between the us. From this, we
        // fallback to an ORect <-> Rect test.
        Oriented< Rect<T,N> > b1_in_b0 = b1 / xf;
        
        return b1_in_b0.intersects(shape);
    }
    
    /// Ray-shape intersection.
    Rect<T,1> intersect(const Ray<T,N>& r) const {
        return shape.intersect(r / xf);
    }
    
    /// Compute the volume of this transformed Rect.
    T volume() const {
        SimpleMatrix<T,N,N> mx;
        for (index_t i = 0; i < N; ++i) {
            // mult this basis vector by the length of the rect along that axis
            T s = shape.hi[i] - shape.lo[i];
            for (index_t j = 0; j < N; ++j) {
                mx(j,i) = s * xf.mat(j,i);
            }
        }
        return det_destructive(mx.data_begin(), N);
    }
    
}; // class Oriented<Rect>


/** @addtogroup shape
 *  @{
 */


/****************************
 * Operators                *
 ****************************/


/**
 * @brief Transform the shape `s` by wrapping it in an `Oriented` class.
 * @related AffineTransform
 * @related Oriented
 */
template <typename Shape>
inline Oriented<Shape> operator*(
        const AffineTransform<typename Shape::elem_t, Shape::N>& xf,
        const Shape& s)
{
    return Oriented<Shape>(s, xf);
}


/**
 * @brief Transform the oriented shape `s` by `xf`.
 * @related AffineTransform
 * @related Oriented
 */
template <typename Shape>
inline Oriented<Shape> operator*(
        const AffineTransform<typename Shape::elem_t, Shape::N>& xf,
        const Oriented<Shape>& s)
{
    return Oriented<Shape>(s.shape, xf * s.xf);
}


/**
 * @brief In-place transform the oriented shape `s` by `xf`.
 * @related AffineTransform
 * @related Oriented
 */
template <typename Shape>
inline Oriented<Shape>& operator*=(
        Oriented<Shape>& s,
        const AffineTransform<typename Shape::elem_t, Shape::N>& xf)
{
    s.xf *= xf;
    return s;
}


/**
 * @brief Transform the shape `s` by the inverse of `xf`.
 * @related AffineTransform
 * @related Oriented
 */
template <typename Shape>
inline Oriented<Shape> operator/(
        const Shape& s,
        const AffineTransform<typename Shape::elem_t, Shape::N>& xf)
{
    return Oriented<Shape>(s, xf.inverse());
}


/**
 * @brief Transform the oriented shape `s` by the inverse of `xf`.
 * @related AffineTransform
 * @related Oriented
 */
template <typename Shape>
inline Oriented<Shape> operator/(
        const Oriented<Shape>& s,
        const AffineTransform<typename Shape::elem_t, Shape::N>& xf)
{
    return Oriented<Shape>(s.shape, s.xf / xf);
}


/**
 * @brief In-place transform the oriented shape `s` by the inverse of `xf`.
 * @related AffineTransform
 * @related Oriented
 */
template <typename Shape>
inline Oriented<Shape>& operator/=(
        Oriented<Shape>& s,
        const AffineTransform<typename Shape::elem_t, Shape::N>& xf)
{
    s.xf /= xf;
    return s;
}

/** @addtogroup traits
 *  @{
 */

// Oriented shapes inherit concepts
template <typename Shape>
struct implements_shape_concept<Oriented<Shape>, RayIntersectable> : 
    public std::integral_constant<
        bool,
        implements_shape_concept<Shape, RayIntersectable>::value>
{};

template <typename Shape>
struct implements_shape_concept<Oriented<Shape>, Convex> : 
    public std::integral_constant<
        bool,
        implements_shape_concept<Shape, RayIntersectable>::value>
{};

/// @} // addtogroup traits
/// @} // addtogroup shape

template <typename Shape, typename H>
struct Digest<Oriented<Shape>, H> {
    H operator()(const Oriented<Shape> &s) const {
        H nonce = geom::truncated_constant<H>(0xb04c31374e530a4e, 0xae9a04ef51d26625);
        return geom::hash_many<H>(nonce, s.shape, s.xf);
    }
};

} // namespace geom


template <typename Shape>
struct std::hash<geom::Oriented<Shape>> {
    size_t operator()(const geom::Oriented<Shape> &s) const {
        return geom::hash<geom::Oriented<Shape>, size_t>(s);
    }
};
