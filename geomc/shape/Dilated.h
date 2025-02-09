#pragma once

#include "geomc/shape/ShapeTypes.h"
#include <geomc/linalg/Similarity.h>
#include <geomc/shape/Rect.h>
#include <geomc/shape/Sphere.h>

/* todo: intersect(ray) could work like this: intersect the dilated box first. if miss,
 * done. the actual hit, if it exists, is inside the interval where the ray overlaps
 * the box. if the ray hits the un-dilated object, this further constrains
 * where the hit could be. then root-find in the interval using sdf().
 *
 * in fact, could use project() instead of sdf() to root-find the ray hit--
 * project() gives you a plane in front of which no shape points are closer.
 * you know immediately that any hit points are behind that plane; which gives
 * you an interval. if at any moment the intersecton of such an interval with any
 * other interval (box interval; another plane interval) is empty, then ray misses.
 * could iteravely project the estimated ray entry/exit points until some kind of
 * convergence.
 *
 * could try to interpolate the hit point— get a point on the inside/outside
 * and interpolate hit `s` value? (would second-order interp help...?)
 *
 * termination conditions:
 * - ds / s?
 * - |surf pt - ray pt|?
 * - ray.v dot n?
 * - |dn|?
 * > want it to work robustly for planar sufraces (i.e. not iterate)
 * 
 *  - another thought was to move the ray origin in the direction toward the surface,
 *    but really you need to move the ray in the direction from the ray axis to the shape.
 *    that's probably not feasible to do without searching, given the methods available.
 *  - second issue is that this wouldn't thicken the shape as we expect.
 */

namespace geom {

/**
 * @ingroup shape
 * @brief A wrapper shape which dilates the extents of another Shape.
 *
 * The shape is the one that would result by placing a sphere of radius `dilation`
 * at every point on the surface of the un-dilated shape.
 *
 * It is not valid to construct a shape with a negative dilation.
 * 
 * Dilated shapes can be constructed with the convenience function `dilate()`,
 * which has overrides for certain shapes that are innately dilatable.
 */
template <typename Shape>
class Dilated: 
    public Convex     <typename Shape::elem_t, Shape::N, Dilated<Shape>>,
    public Projectable<typename Shape::elem_t, Shape::N, Dilated<Shape>>
{
public:
    /// The coordinate type of this Shape
    typedef typename Shape::elem_t T;
    /// The dimension of this Shape
    static constexpr size_t N = Shape::N;
    
    /// Un-dilated shape.
    Shape shape;
    /// Dilation distance; a nonnegative number.
    T dilation;
    
    /// Wrap a default-constructed Shape with zero dilation.
    Dilated():dilation((T)0) {}
    
    /// Wrap `s` with zero dilation.
    Dilated(const Shape& s):shape(s) {}
    
    /// Dilate `s` by the amount `dilation`.
    Dilated(const Shape& s, T dilation):
        shape(s),
        dilation(std::max<T>(dilation, 0)) {}
        
    static constexpr bool admits_cusps() { return false; }
    
    bool operator==(const Dilated& other) const {
        return shape == other.shape && dilation == other.dilation;
    }
    
    /// Point containment test.
    bool contains(Vec<T,N> p) const requires RegionObject<Shape> {
        return sdf(p) < 0;
    }
    
    /// Signed distance function.
    T sdf(Vec<T,N> p) const requires SdfObject<Shape> {
        return shape.sdf(p) - dilation;
    }
    
    Vec<T,N> normal(Vec<T,N> p) const requires ProjectableObject<Shape> {
        // we are dilating the shape along the normal, so the normal is unchanged
        return shape.normal(p);
    }
    
    Vec<T,N> convex_support(Vec<T,N> d) const requires ConvexObject<Shape> {
        Vec<T,N> p = shape.convex_support(d);
        return p + d.unit() * dilation;
    }
    
    Rect<T,N> bounds() const requires BoundedObject<Shape> {
        return shape.bounds().dilated(dilation);
    }
    
    /// Orthogonally project `p` to the surface of this shape.
    Vec<T,N> project(Vec<T,N> p) const
            requires ProjectableObject<Shape> and ConvexObject<Shape>
    {
        Vec<T,N> p_proj = shape.project(p);
        Vec<T,N> dp = shape.normal(p);
        if (shape.contains(p)) {
            // dilated projection is further away from the point
            return p_proj + dilation * dp;
        } else {
            // dilated projection is more toward point
            return p_proj - dilation * dp;
        }
    }
    
    // todo: trace()
};


/** @addtogroup shape
 *  @{
 */

/**
 * @brief Dilate the shape `s` by the amount `dilation`.
 * @related Dilated
 */
template <typename Shape>
inline Dilated<Shape> dilate(const Shape& s, typename Shape::elem_t dilation) {
    return Dilated<Shape>(s, std::max(dilation, 0));
}

/**
 * @brief Dilate the shape `s` by the amount `dilation`.
 * @related Dilated
 */
template <typename Shape>
inline Dilated<Shape> dilate(const Dilated<Shape>& s, typename Shape::elem_t dilation) {
    return Dilated<Shape>(s.shape, std::max(s.dilation + dilation, 0));
}

/**
 * @brief Dilate the Plane `p` by the amount `dilation`.
 *
 * (This just moves the plane `dilation` units along its normal).
 * @related Plane
 * @related Dilated
 */
template <typename T, index_t N>
inline Plane<T,N> dilate(Plane<T,N> p, T dilation) {
    p.d += dilation;
    return p;
}

/**
 * @brief Create a spherical shell around `center`, with inner and outer radii `r0` and `r1`.
 *
 * It is not required that r0 < r1; the shell will always be constructed with positive thickness.
 * @related Dilated
 */
template <typename T, index_t N>
inline Dilated<Hollow<Sphere<T,N>>> shell(Vec<T,N> center, T r0, T r1) {
    r0 = std::max<T>(r0, 0);
    r1 = std::max<T>(r1, 0);
    return {
        Sphere<T,N>{center, (r0 + r1) / 2},
        std::abs(r1 - r0) / 2
    };
}

/**
 * @brief Create a spherical shell with wall thickness `thickness` around the Sphere `s`.
 *
 * The wall will be centered on the Sphere's surface.
 * @related Dilated
 */
template <typename T, index_t N>
inline Dilated<Hollow<Sphere<T,N>>> shell(const Sphere<T,N>& s, T thickness) {
    T t = thickness / 2;
    return shell(s.center, s.r - t, s.r + t);
}

/**
 * @brief Dilate the Sphere `s` by the amount `dilation`.
 *
 * (This just adds the amount `dilation` to the Sphere's radius).
 * @related Sphere
 * @related Dilated
 */
template <typename T, index_t N>
inline Sphere<T,N> dilate(const Sphere<T,N>& s, T dilation) {
    return Sphere<T,N>(s.center, std::max(s.r + dilation, 0));
}

/**
 * @brief Produce a rounded rectangle with corner radius `radius`.
 * 
 * The rounded rectangle will have the extents of `rect`, but with corners
 * rounded by `radius`. The radius will be clamped to the maximum possible
 * value that will not cause the corners to overlap.
 */
template <typename T, index_t N>
inline Dilated<Rect<T,N>> roundrect(const Rect<T,N>& rect, T radius) {
    radius = std::min(radius, rect.dimensions().min() / 2);
    return {rect.dilated(-radius), radius};
}

/// @brief Transform a dilated shape by a similarity transform.
/// @related Dilated
/// @related Similarity
template <typename T, index_t N, typename Shape>
requires (N == Shape::N and std::same_as<T, typename Shape::elem_t> and
    requires (Shape s, Similarity<T,N> xf) { {s * xf} -> std::same_as<Shape>; }
)
Dilated<Shape> operator*(const Similarity<T,N>& xf, const Dilated<Shape>& s) {
    return Dilated<Shape>(xf * s.shape, s.dilation * xf.sx);
}

/// @brief Inverse-transform a dilated shape by a similarity transform.
/// @related Dilated
/// @related Similarity
template <typename T, index_t N, typename Shape>
requires (N == Shape::N and std::same_as<T, typename Shape::elem_t> and
    requires (Shape s, Similarity<T,N> xf) { {s / xf} -> std::same_as<Shape>; }
)
Dilated<Shape> operator/(const Dilated<Shape>& s, const Similarity<T,N>& xf) {
    return Dilated<Shape>(xf / s.shape, s.dilation / xf.sx);
}

/// @brief Transform a dilated shape by an isometry.
/// @related Dilated
/// @related Isometry
template <typename T, index_t N, typename Shape>
requires (N == Shape::N and std::same_as<T, typename Shape::elem_t> and
    requires (Shape s, Isometry<T,N> xf) { {s * xf} -> std::same_as<Shape>; }
)
Dilated<Shape> operator*(const Isometry<T,N>& xf, const Dilated<Shape>& s) {
    return Dilated<Shape>(xf * s.shape, s.dilation);
}

/// @brief Inverse-transform a dilated shape by an isometry.
/// @related Dilated
/// @related Isometry
template <typename T, index_t N, typename Shape>
requires (N == Shape::N and std::same_as<T, typename Shape::elem_t> and
    requires (Shape s, Isometry<T,N> xf) { {s / xf} -> std::same_as<Shape>; }
)
Dilated<Shape> operator/(const Dilated<Shape>& s, const Isometry<T,N>& xf) {
    return Dilated<Shape>(xf / s.shape, s.dilation);
}

/// @} // group shape

template <typename Shape, typename H>
struct Digest<Dilated<Shape>, H> {
    H operator()(const Dilated<Shape> &s) const {
        H nonce = geom::truncated_constant<H>(0x7f5e38d0b18f1dbf, 0x105e0a188e1bc0ef);
        return geom::hash_many<H>(nonce, s.shape, s.dilation);
    }
};

} // namespace geom


template <typename Shape>
struct std::hash<geom::Dilated<Shape>> {
    size_t operator()(const geom::Dilated<Shape> &s) const {
        return geom::hash<geom::Dilated<Shape>, size_t>(s);
    }
};
