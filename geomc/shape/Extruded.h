#pragma once

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
 *     auto cylinder = Extruded<Sphere<double,2>>(Sphere<double,2>(), Rect<double,1>(0,1));
 *     // or equivalently, and more succinctly:
 *     auto cylinder2 = extrude(Circle<double>(), 0, 1);
 */
template <typename Shape>
class Extruded: public Dimensional<typename Shape::elem_t, Shape::N + 1> {
public:
    using typename Dimensional<typename Shape::elem_t, Shape::N + 1>::elem_t;
    using Dimensional<typename Shape::elem_t, Shape::N + 1>::N;
    using T = elem_t;
    
    /// Cross-section of this extrusion at the h = 0 plane.
    Shape base;
    /// Height range of this extrusion.
    Rect<T,1> height;
    
    
    /// Construct a new extrusion with unit height.
    Extruded():height(0,1) {}
    
    /// Construct a new extrusion with cross section `base` and unit height.
    explicit Extruded(const Shape& base):
        base(base),
        height(0,1) {}
    
    /// Construct an extrusion with cross section `base` and height range `height`.
    Extruded(const Shape& base, const Rect<T,1>& height):
            base(base), 
            height(height) {}
    
    /** @brief Construct an extrusion with cross section `base`, and height 
     * ranging between `h0` and `h1`.
     */     
    Extruded(const Shape& base, T h0, T h1):
        base(base),
        height(std::min(h0, h1), std::max(h0, h1)) {}
        
    static constexpr bool admits_cusps() { return Shape::admits_cusps(); }
    
    bool operator==(const Extruded<Shape>& other) const {
        return base == other.base && height == other.height;
    }
    
    bool contains(Vec<T,N> p) const requires RegionObject<Shape> {
        if (not height.contains(p[N-1])) return false;
        return base.contains(p.template resized<N-1>());
    }
    
    T sdf(Vec<T,N> p) const requires SdfObject<Shape> {
        Vec<T,N-1> p_proj = (Vec<T,N-1>) p.template resized<N-1>();
        T h = p[N-1];
        // distance to the base shape within the base shape's plane
        T sdf_base = base.sdf(p_proj);
        // sdf of the slab that contains the height range of the shape
        T sdf_h_slab = std::max(h - height.hi, height.lo - h);
        if (height.contains(h)) {
            // intersect base extrusion with the height slab
            // because we are inside the shape, the sdf will be negative,
            // so max() will give us a smaller absolute value— the closer surface.
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
    
    Vec<T,N> normal(Vec<T,N> p) const requires ProjectableObject<Shape> {
        Vec<T,N-1> p_base = (Vec<T,N-1>) p.template resized<N-1>();
        T h = p[N-1];
        bool inside_base = base.contains(p_base);
        if (height.contains(h)) {
            // p is somewhere between the caps
            if (inside_base) {
                // we need to know if the wall is closest, or a cap is closest
                T cap_h     = height.project(h);
                T cap_dist  = std::abs(h - cap_h);
                T wall_dist = std::abs(base.sdf(p_base));
                if (cap_dist < wall_dist) {
                    // the cap is nearest. use its normal
                    return {Vec<T,N-1>{}, (cap_h >= height.hi) ? 1 : -1};
                }
            }
            // the wall is nearest, use its normal
            return Vec<T,N>(base.normal(p_base), 0);
        } else {
            // p beyond the top or bottom cap
            if (inside_base) {
                // p directly above or below a cap
                return {Vec<T,N-1>{}, (h >= height.hi) ? 1 : -1};
            } else {
                // p is outside the base shape
                Vec<T,N> cap_pt {
                    base.project(p_base),
                    height.clip(h)
                };
                return (p - cap_pt).unit();
            }
        }
    }
    
    /**
     * Return the point `p` orthogonally projected onto the surface of the shape.
     */
    Vec<T,N> project(Vec<T,N> p) const requires ProjectableObject<Shape> {
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
    
    Vec<T,N> clip(Vec<T,N> p) const requires ProjectableObject<Shape> {
        return contains(p) ? p : project(p);
    }
    
    Vec<T,N> convex_support(Vec<T,N> d) const requires ConvexObject<Shape> {
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
    
    Rect<T,N> bounds() const requires BoundedObject<Shape> {
        return base.bounds() * height;
    }
    
    template <ConvexObject S>
    requires ConvexObject<Shape> and (S::N == N) and std::same_as<T, typename S::elem_t>
    bool intersects(const S& other) const {
        return geom::intersects(
            as_any_convex(*this),
            as_any_convex(other)
        );
    }
    
    /// Measure the interior (volume) of the extrusion.
    T measure_interior() const requires InteriorMeasurableObject<Shape> {
        return base.measure_interior() * height.dimensions();
    }
    
    /// Measure the boundary (surface area) of the extrusion.
    T measure_boundary() const requires BoundaryMeasurableObject<Shape> {
        T cap_area = 0;
        // it is possible to extrude a boundary shape to make a "tube";
        // in that case it will not have caps. only count the cap area if
        // the shape has an interior.
        if constexpr (InteriorMeasurableObject<Shape>) {
            cap_area = 2 * base.measure_interior();
        }
        return (
            // walls
            base.measure_boundary() * height.dimensions() +
            // caps
            cap_area
        );
    }
    
    /// Ray/shape intersection.
    Rect<T,1> intersect(const Ray<T,N>& r) const requires RayIntersectableObject<Shape> {
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
            interval &= Rect<T,1>::from_corners(
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
 * wrapping `s` in the `Extruded` template.
 *
 * @related Extruded
 */
template <typename Shape>
inline Extruded<Shape> extrude(
        const    Shape& s,
        typename Shape::elem_t h0,
        typename Shape::elem_t h1)
{
    return Extruded<Shape>(s, h0, h1);
}

/// @} // addtogroup shape


template <typename Shape, typename H>
struct Digest<Extruded<Shape>, H> {
    H operator()(const Extruded<Shape> &s) const {
        H nonce = geom::truncated_constant<H>(0x28211b7d8ba5f09b, 0x6b7b5d58d68273c7);
        return geom::hash_many<H>(nonce, s.base, s.height);
    }
};

#ifdef GEOMC_USE_STREAMS

template <typename Shape>
std::ostream& operator<<(std::ostream& os, const Extruded<Shape>& e) {
    os << "Extruded(" << e.shape << ", " << e.height << ")";
    return os;
}

#endif

} // namespace geom


template <typename Shape>
struct std::hash<geom::Extruded<Shape>> {
    size_t operator()(const geom::Extruded<Shape> &s) const {
        return geom::hash<geom::Extruded<Shape>, size_t>(s);
    }
};
