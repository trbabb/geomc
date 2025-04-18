#pragma once

#include <geomc/linalg/Similarity.h>
#include <geomc/shape/Rect.h>

// todo: add intersection for SAT-tested bboxes

namespace geom {

/**
 * @ingroup shape
 * @brief A shape transformed by a similarity transform (translation, rotation, and scale).
 * 
 * Similar transfrorms guarantee that shapes and relative distances are preserved.
 * For arbitrarily-transformed shapes, see Transformed.
 */
template <typename Shape>
class Similar: public Dimensional<typename Shape::elem_t, Shape::N> {
public:
    using elem_t = typename Shape::elem_t;
    static constexpr index_t N = Shape::N;
    
    using T = elem_t;
    
    /// Un-transformed shape
    Shape shape;
    Similarity<T,N> xf;
    
    /// Wrap a default-constructed shape with the identity transform.
    Similar() {}
    
    /// Wrap a shape with a similarity transform.
    Similar(const Shape& shape, const Similarity<T,N>& xf):
        shape(shape),
        xf(xf) {}
    
    /// Wrap a shape with an identity transform.
    Similar(const Shape& shape):
        shape(shape) {}
        
    static constexpr bool admits_cusps() { return Shape::admits_cusps(); }
    
    bool operator==(const Similar& other) const {
        return shape == other.shape && xf == other.xf;
    }
    
    /// Shape-point intersection test.
    bool contains(typename Shape::point_t p) const requires RegionObject<Shape> {
        return shape.contains(p / xf);
    }
    
    /// Convert this shape to a Transformed shape.
    operator Transformed<Shape>() const {
        return Transformed<Shape>(shape, xf);
    }
    
    template <ConvexObject S>
    bool intersects(const S& other) const {
        return geom::intersects(
            as_any_convex(*this),
            as_any_convex(other)
        );
    }
    
    /// @brief Intersecion with another similar shape.
    /// Available if our shape can intersect the other shape's base shape.
    template <typename S>
    bool intersects(const Similar<S>& s) const 
        requires requires (const Similar<Shape> s, const S b) { s.intersects(b); }
    {
        // if we can intersect with an un-transformed shape,
        // then we can also intersect with a transformed one by
        // un-transforming ourselves first.
        Similar<Shape> s1 = *this / s.xf;
        return s1.intersects(s.shape);
    }
    
    /// @brief Intersection with another shape.
    /// Available if the other shape can be transformed into our space
    /// while preserving its type.
    template <typename S>
    bool intersects(const S& s) const 
        requires requires (Similarity<T,N> xf, const Shape s, const S b) { 
            s.intersects(b);
            { b / xf } -> std::same_as<S>;
        }
    {
        // if we can transform the other shape into base-shape space,
        // and our base shape can intersect the other shape,
        // then we can intersect with the other shape.
        return shape.intersects(s / xf);
    }
    
    /// Convex support function. Return the point on the surface of the shape
    /// which is farthest in the direction `d`.
    Vec<T,N> convex_support(Vec<T,N> d) const requires ConvexObject<Shape> {
        return xf * shape.convex_support(xf.apply_inverse_direction(d));
    }
    
    /// Compute the axis-aligned bounding box of the shape.
    Rect<T,N> bounds() const requires BoundedObject<Shape> {
        Rect<T,N> r;
        for (index_t axis = 0; axis < N; ++axis) {
            // test along a cardinal axis in shape-space
            Vec<T,N> d;
            d[axis] = 1;
            d = xf.apply_inverse_direction(d);
            Vec<T,N> p_lo = xf * shape.convex_support( d);
            Vec<T,N> p_hi = xf * shape.convex_support(-d);
            r.lo[axis] = p_lo[axis];
            r.hi[axis] = p_hi[axis];
        }
        return r;
    }
    
    Vec<T,N> corner(index_t i) const requires std::same_as<Shape,Rect<T,N>> {
        return xf * shape.corner(i);
    }
    
    Similar<Rect<T,N>> transformed_bounds() const requires BoundedObject<Shape> {
        return Similar<Rect<T,N>>(bounds(), xf);
    }
    
    /// Ray-shape intersection.
    Rect<T,1> intersect(const Ray<T,N>& ray) const requires RayIntersectableObject<Shape> {
        return xf * shape.intersect(ray / xf);
    }
    
    /// Signed distance function.
    T sdf(Vec<T,N> p) const requires SdfObject<Shape> {
        return xf.sx * shape.sdf(p / xf);
    }
    
    /// Direction away from the surface of the shape at point `p`.
    Vec<T,N> normal(Vec<T,N> p) const requires ProjectableObject<Shape> {
        return xf.rx * shape.normal(p / xf);
    }
    
    /// Orthogonally project `p` to the surface of this shape.
    Vec<T,N> project(Vec<T,N> p) const requires ProjectableObject<Shape> {
        return xf * shape.project(p / xf);
    }
    
    /// Nearest point on the interior of the shape.
    Vec<T,N> clip(Vec<T,N> p) const requires ProjectableObject<Shape> {
        return contains(p) ? p : project(p);
    }
    
    /// Measure the interior (volume) of the shape.
    T measure_interior() const requires InteriorMeasurableObject<Shape> {
        return xf.sx * shape.measure_interior();
    }
    
    /// Measure the boundary (surface area) of the shape.
    T measure_boundary() const requires BoundaryMeasurableObject<Shape> {
        return xf.sx * shape.measure_boundary();
    }
    
};

/// @addtogroup shape
/// @{

/// @brief Transform the shape `shape` by wrapping it with a Similarity transform.
/// @related Similar
template <typename Shape>
inline Similar<Shape> operator*(
        const Similarity<typename Shape::elem_t, Shape::N>& xf,
        const Shape& shape)
{
    return Similar<Shape>(shape, xf);
}

/// @brief Transform the shape `s` by `xf`.
/// @related Similar
template <typename Shape>
inline Similar<Shape> operator*(
    const Similarity<typename Shape::elem_t, Shape::N>& xf,
    const Similar<Shape>& s)
{
    return Similar<Shape>(s.shape, xf * s.xf);
}

/// @brief In-place transform the shape `s` by `xf`.
/// @related Similar
template <typename Shape>
inline Similar<Shape>& operator*=(
    Similar<Shape>& s,
    const Similarity<typename Shape::elem_t, Shape::N>& xf)
{
    s.xf = xf * s.xf;
    return s;
}

/// @brief Transform the shape `s` by the inverse of `xf`.
/// @related Similar
template <typename Shape>
inline Similar<Shape> operator/(
    const Shape& s,
    const Similarity<typename Shape::elem_t, Shape::N>& xf)
{
    return Similar<Shape>(s, xf.inv());
}

/// @brief Transform the shape `s` by the inverse of `xf`.
/// @related Similar
template <typename Shape>
inline Similar<Shape> operator/(
    const Similar<Shape>& s,
    const Similarity<typename Shape::elem_t, Shape::N>& xf)
{
    return Similar<Shape>(s.shape, s.xf / xf);
}

/// @brief In-place transform the shape `s` by the inverse of `xf`.
/// @related Similar
template <typename Shape>
inline Similar<Shape>& operator/=(
    Similar<Shape>& s,
    const Similarity<typename Shape::elem_t, Shape::N>& xf)
{
    s.xf = s.xf / xf;
    return s;
}

/// @} // group shape


template <typename Shape, typename H>
struct Digest<Similar<Shape>, H> {
    H operator()(const Similar<Shape>& s, H& h) const {
        H nonce = geom::truncated_constant<H>(0x742da870d5a73569, 0xbbec1fb3638d6150);
        return geom::hash_many<H>(nonce, s.shape, s.xf);
    }
};

#ifdef GEOMC_USE_STREAMS

template <typename Shape>
std::ostream& operator<<(std::ostream& os, const Similar<Shape>& s) {
    os << "Similar(" << s.shape << ", " << s.xf << ")";
    return os;
}

#endif

} // namespace geom


template <typename Shape>
struct std::hash<geom::Similar<Shape>> {
    size_t operator()(const geom::Similar<Shape> &s) const {
        return geom::hash<geom::Similar<Shape>, size_t>(s);
    }
};
