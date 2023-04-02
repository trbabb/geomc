#pragma once
/* 
 * File:   Frustum.h
 * Author: tbabb
 *
 * Created on November 8, 2014, 11:09 PM
 */

#include <climits>
#include <geomc/shape/Rect.h>
#include <geomc/linalg/AffineTransform.h>


namespace geom {

/**
 * @ingroup shape
 * @brief An axis-aligned extrusion of an arbitrary N-1 dimensional Convex shape.
 * 
 * The first N-1 dimensions have cross-sections which are `Shape`s
 * (e.g., `Rect`, `Sphere`, `Frustum`...). The shape is lofted "vertically" into 
 * the last dimension.
 *
 * Example:
 * 
 *     // initialize a cylinder of length and radius 1:
 *     auto cylinder = Extrusion<Sphere<double,2>>(Sphere<double,2>(), Rect<double,1>(0,1));
 *     // or equivalently, and more succinctly:
 *     auto cylinder2 = extrude(Circle<double>(), 0, 1);
 */
template <typename Shape>
class Extrusion :
    public Convex          <typename Shape::elem_t, Shape::N+1, Extrusion<Shape>>,
    public RayIntersectable<typename Shape::elem_t, Shape::N+1, Extrusion<Shape>>,
    public Projectable     <typename Shape::elem_t, Shape::N+1, Extrusion<Shape>>
{
public:
    /// The coordinate type of this Shape
    typedef typename Shape::elem_t T;
    /// The dimension of this Shape
    static constexpr size_t N = Shape::N + 1;
    
    /// Cross-section of this extrusion at the h = 0 plane.
    Shape base;
    /// Height range of this extrusion.
    Rect<T,1> height;
    
    
    /// Construct a new extrusion with unit height.
    Extrusion():height(0,1) {}
    
    /// Construct a new extrusion with cross section `base` and unit height.
    explicit Extrusion(const Shape& base):
        base(base),
        height(0,1) {}
    
    /// Construct an extrusion with cross section `base` and height range `height`.
    Extrusion(const Shape& base, const Rect<T,1>& height):
            base(base), 
            height(height) {}
    
    /** @brief Construct an extrusion with cross section `base`, and height 
     * ranging between `h0` and `h1`.
     */     
    Extrusion(const Shape& base, T h0, T h1):
        base(base),
        height(std::min(h0, h1), std::max(h0, h1)) {}
    
    
    bool contains(Vec<T,N> p) const {
        if (not height.contains(p[N-1])) return false;
        return base.contains(p.template resized<N-1>());
    }
    
    inline T sdf(Vec<T,N> p) const {
        Vec<T,N-1> p_proj = (Vec<T,N-1>) p.template resized<N-1>();
        T h = p[N-1];
        // distance to the base shape within the base shape's plane
        T sdf_base = base.sdf(p_proj);
        // sdf of the slab that contains the height range of the shape
        T sdf_h_slab = std::max(h - height.hi, height.lo - h);
        if (height.contains((T)h)) {
            // intersect base extrusion with the height slab
            return std::max(sdf_base, sdf_h_slab);
        } else {
            // nearest point must be a cap
            return (sdf_base < 0)
                // point projects orthogonally to the cap
                ? sdf_h_slab
                // point is "away" from the edge of the cap, above/below it
                : std::sqrt(sdf_base * sdf_base + sdf_h_slab * sdf_h_slab);
        }
    }
    
    /**
     * Return the point `p` orthogonally projected onto the surface of the shape.
     */
    inline Vec<T,N> project(Vec<T,N> p) const {
        const T h     = p[N-1];
        const T h_cap = std::abs(h - height.lo) < std::abs(h - height.hi) 
            ? height.lo
            : height.hi;
        // point in the plane of the base shape:
        const Vec<T,N-1> p_base = p.template resized<N-1>();
        // point projected to the boundary of the base shape:
        const Vec<T,N-1> p_proj = base.project(p_base);
        // nearest point on the wall:
        Vec<T,N> wall_pt = Vec<T,N>(p_proj, height.clip(h));
        // point on the nearer cap:
        Vec<T,N> cap_pt  = base.contains(p_base)
            ? Vec<T,N>(p_base, h_cap)
            : Vec<T,N>(p_proj, h_cap);
        // pick the closer of the cap or the wall:
        return p.dist2(wall_pt) < p.dist2(cap_pt) ? wall_pt : cap_pt;
    }
    
    Vec<T,N> convex_support(Vec<T,N> d) const {
        typedef PointType<T,N-1> Pt;
        
        Vec<T,N-1> d_ = d.template resized<N-1>();
        
        if (d_.is_zero()) {
            // support is normal to the endcap
            // pick an arbitrary point on the endcap
            d_[0] = 1;
        }
        
        Vec<T,N-1> p0 = base.convex_support(d_);
        
        // top or bottom face?
        return Vec<T,N>(p0, (d[N-1] > 0) ? height.hi : height.lo);
    }
    
    Rect<T,N> bounds() const {
        return base.bounds() * height;
    }
    
    /// Ray/shape intersection.
    Rect<T,1> intersect(const Ray<T,N>& r) const {
        Ray<T,N-1> r_base = r.template resized<N-1>();
        Rect<T,1> interval;
        if (r_base.direction.is_zero()) {
            // ray is pointing directly along extrusion axis
            if (base.contains(r_base.origin)) {
                interval = Rect<T,1>::full;
            }
        } else {
            interval = base.intersect(r_base);
        }
        T h  = r.origin[N-1];
        T dH = r.direction[N-1];
        if (dH == 0) {
            // ray is parallel to the base plane
            return height.contains(h) ? interval : Rect<T,1>::empty;
        } else {
            // trace the height slab and intersect its interval
            interval &= Rect<T,1>::spanning_corners(
                (height.lo - h) / dH,
                (height.hi - h) / dH);
            return interval;
        }
    }


}; // class extrusion


/** 
 * @addtogroup shape
 * @{
 */

/**
 * @brief Convenience function to extrude the shape `s` between heights `h0` and `h1` by
 * wrapping `s` in the `Extrusion` template.
 *
 * @related Extrusion
 */
template <typename Shape>
inline Extrusion<Shape> extrude(
        const    Shape& s,
        typename Shape::elem_t h0,
        typename Shape::elem_t h1)
{
    return Extrusion<Shape>(s, h0, h1);
}

/** 
 * @addtogroup traits 
 * @{
 */

// Extruded shapes inherit concepts
template <typename Shape>
struct implements_shape_concept<Extrusion<Shape>, Projectable> : 
    public std::integral_constant<
        bool,
        implements_shape_concept<Shape, Projectable>::value>
{};

template <typename Shape>
struct implements_shape_concept<Extrusion<Shape>, RayIntersectable> : 
    public std::integral_constant<
        bool,
        implements_shape_concept<Shape, RayIntersectable>::value>
{};

template <typename Shape>
struct implements_shape_concept<Extrusion<Shape>, Convex> : 
    public std::integral_constant<
        bool,
        implements_shape_concept<Shape, RayIntersectable>::value>
{};

template <typename Shape>
struct implements_shape_concept<Extrusion<Shape>, SdfEvaluable> : 
    public std::integral_constant<
        bool,
        implements_shape_concept<Shape, SdfEvaluable>::value>
{};

/// @} // addtogroup traits
/// @} // addtogroup shape

} // namespace geom
