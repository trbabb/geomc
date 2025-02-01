#pragma once

#include <geomc/linalg/Vec.h>
#include <geomc/linalg/Similarity.h>
#include <geomc/shape/Shape.h>

namespace geom {

/**
 * @ingroup shape
 * @brief Selects the boundary of a shape.
 * 
 * Hollow shapes are infinitely thin. They are the boundary of a shape, but do not
 * contain any points. Hollow shapes can report the distance to their surface,
 * the normal vector at a test point, and can orthogonally project points onto
 * their surface.
 *
 * Hollow shapes may be a useful foundation for thick-shelled shapes using
 * the `Dilated` shape wrapper.
 */
template <typename Shape>
class Hollow:
    public Projectable<typename Shape::elem_t, Shape::N, Hollow<Shape>>
    // - hollow shapes are in general not convex
    // - hollow shapes do not have one ray-intersection interval
{
public:
    using elem_t = typename Shape::elem_t;
    static constexpr index_t N = Shape::N;
    
    using T = elem_t;
    
    /// Solid shape.
    Shape shape;
    
    /// Wrap a default-constructed shape with the identity transform.
    Hollow() {}
    
    /// Wrap a shape with a similarity transform.
    Hollow(const Shape& shape):
        shape(shape) {}
    
    /// Shape equality test.
    bool operator==(const Hollow& other) const {
        return shape == other.shape;
    }
    
    /**
     * @brief Shape-point intersection test. Always returns `false`, because 
     * hollow shapes are infinitely thin.
     */
    bool contains(typename Shape::point_t p) const {
        return false; // hollow shapes are infinitely thin
    }
    
    /// Project a point to the surface of this shape.
    Vec<T,N> project(Vec<T,N> p) const {
        return shape.project(p);
    }
    
    /**
     * @brief Nearest point on the interior of the shape to `p`; the same as
     * `project(p)`, since hollow shapes do not have an interior.
     */
    Vec<T,N> clip(Vec<T,N> p) const {
        return shape.project(p);
    }
    
    /// Normal vector at point `p`.
    Vec<T,N> normal(Vec<T,N> p) const {
        T sign = shape.contains(p) ? -1 : 1;
        return sign * shape.normal(p);
    }
    
    /// Signed distance function.
    T sdf(Vec<T,N> p) const {
        return std::abs(shape.sdf(p));
    }
    
    /// Compute the axis-aligned bounding box of the shape.
    Rect<T,N> bounds() const {
        return shape.bounds();
    }
    
};

/// @ingroup shape
/// @{

/**
 * @brief Transform a hollow shape.
 * @related Hollow
 */
template <typename Shape, Transform<typename Shape::elem_t, Shape::N> Xf>
requires Transformable<Shape, Xf>
inline Hollow<Shape> operator*(const Xf& xf, const Hollow<Shape>& s) {
    return Hollow<Shape>(xf * s.shape);
}

/**
 * @brief Inverse transform a hollow shape.
 * @related Hollow
 */
template <typename Shape, Transform<typename Shape::elem_t, Shape::N> Xf>
requires Transformable<Shape, Xf>
inline Hollow<Shape> operator/(const Hollow<Shape>& s, const Xf& xf) {
    return Hollow<Shape>(xf / s.shape);
}

/**
 * @brief In-place transform a hollow shape.
 * @related Hollow
 */
template <typename Shape, Transform<typename Shape::elem_t, Shape::N> Xf>
requires Transformable<Shape, Xf>
inline Hollow<Shape>& operator*=(Hollow<Shape>& s, const Xf& xf) {
    s.shape *= xf;
    return s;
}

/**
 * @brief In-place inverse transform a hollow shape.
 * @related Hollow
 */
template <typename Shape, Transform<typename Shape::elem_t, Shape::N> Xf>
requires Transformable<Shape, Xf>
inline Hollow<Shape>& operator/=(Hollow<Shape>& s, const Xf& xf) {
    s.shape /= xf;
    return s;
}

/// @}

template <typename Shape, typename H>
struct Digest<Hollow<Shape>, H> {
    H operator()(const geom::Hollow<Shape> &s) const {
        H nonce = geom::truncated_constant<H>(0xecd1792decf07dcb, 0xdd9a7673d74739b2);
        return geom::hash_many(nonce, s.shape);
    }
};

}  // namespace geom

template <typename Shape>
struct std::hash<geom::Hollow<Shape>> {
    size_t operator()(const geom::Hollow<Shape> &s) const {
        return geom::hash<geom::Hollow<Shape>, size_t>(s);
    }
};
