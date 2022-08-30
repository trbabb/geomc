#ifndef _DILATED_SHAPE_H
#define _DILATED_SHAPE_H

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
    Dilated(const Shape& s, T dilation):shape(s), dilation(dilation) {}
    
    /// Point containment test.
    bool contains(Vec<T,N> p) const {
        return sdf(p) < 0;
    }
    
    /// Signed distance function.
    inline T sdf(Vec<T,N> p) const {
        return shape.sdf(p) - dilation;
    }
    
    Vec<T,N> convex_support(Vec<T,N> d) const {
        Vec<T,N> p = shape.convex_support(d);
        return p + d.unit() * dilation;
    }
    
    Rect<T,N> bounds() const {
        return shape.bounds().dilated(dilation);
    }
    
    /// Orthogonally project `p` to the surface of this shape.
    Vec<T,N> project(Vec<T,N> p) const {
        Vec<T,N> p_proj = shape.project(p);
        Vec<T,N> dp = (p_proj - p).unit();
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
Dilated<Shape> dilate(const Shape& s, typename Shape::elem_t dilation) {
    return Dilated<Shape>(s, std::max(dilation, 0));
}

/**
 * @brief Dilate the shape `s` by the amount `dilation`.
 * @related Dilated
 */
template <typename Shape>
Dilated<Shape> dilate(const Dilated<Shape>& s, typename Shape::elem_t dilation) {
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
Plane<T,N> dilate(Plane<T,N> p, T dilation) {
    p.d += dilation;
    return p;
}

/**
 * @brief Dilate the Sphere `s` by the amount `dilation`.
 *
 * (This just adds the amount `dilation` to the Sphere's radius).
 * @related Sphere
 * @related Dilated
 */
template <typename T, index_t N>
Sphere<T,N> dilate(const Sphere<T,N>& s, T dilation) {
    return Sphere<T,N>(s.center, std::max(s.r + dilation, 0));
}

/** @addtogroup traits
 *  @{
 */

// Dilated shapes inherit concepts
template <typename Shape>
struct implements_shape_concept<Dilated<Shape>, Projectable> : 
    public std::integral_constant<
        bool,
        implements_shape_concept<Shape, Projectable>::value>
{};

template <typename Shape>
struct implements_shape_concept<Dilated<Shape>, Convex> : 
    public std::integral_constant<
        bool,
        implements_shape_concept<Shape, RayIntersectable>::value>
{};

/// @}  // addtogroup traits
/// @}  // addtogroup shape

} // namespace geom

#endif /* _DILATED_SHAPE_H
