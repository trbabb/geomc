#pragma once

#include <geomc/shape/Rect.h>
#include <geomc/shape/Sphere.h>

namespace geom {

// todo: point projection and sdf is unsolved
//   - can't project point to base plane and recurse
//   - example case: consider a shape far from origin
//     - pt will project to the wall directly "underneath" it, in the limit
//     - it's not just the Z component that's distorted; the direction to the shape
//       surface is also distorted
//     - in general, projection vector field is "stretched" away from origin by this effect
//   - answer: we can perform the projection if we can project to the surface of a 
//     *perspective transformation* of the base shape. that is, put the base shape on a piece
//     of glass, shine a point light at it, and let the shadow fall onto an arbitrarily-oriented
//     plane. We want to project to the shape of the shadow.
//     - the plane we are projecting to is sensitive to the point we are testing. that is, the
//       plane is the one normal to the test point vector (i.e., a vector with its tail at
//       the origin / frustum tip). proof is too small for margin to contain.
//     - sdf() draws a circle around the point which contains none of the shape;
//       this circle becomes a conic section in the perspective-transformed space. the
//       actual closest point may be inside the circle that circumscribes the conic section;
//       we do not know in which direction that point lies from sdf() or project() alone.
//     - in general this cannot be written in terms of the SDF or projection of the base shape
//       without some kind of iterative solution.

// todo: permit 1d cross sectional shapes

/**
 * @ingroup shape
 * @brief An N-dimensional frustum (truncated pyramid) with an arbitrary
 * Convex shape as its base, and its (possibly excluded) point at the origin.
 * 
 * The first N-1 dimensions have cross-sections which are `Shape`s. The last
 * dimension is the "height" of the frustum. The base `Shape` is specified 
 * lying on the `h = 1` plane. (Note that below the origin, the base shape
 * will be flipped due to the change of sign).
 *  
 * If the height range spans `height = 0`, then the height extending below the 
 * origin is excluded from the shape. Use `clipped_height()` to obtain the 
 * height range actually spanned by the shape. (A frustum spanning the origin
 * would not be convex).
 *
 * Oriented, Rect-based Frustums are very commonly needed to represent
 * viewing frustums; for these consider using `ViewFrustum`, which is a
 * templated type alias for `Transformed<Frustum<Rect>>`. 
 */
template <typename Shape>
class Frustum: public Dimensional<typename Shape::elem_t, Shape::N+1> {
public:
    /// The coordinate type of this Shape
    typedef typename Shape::elem_t T;
    /// The dimension of this Shape
    static constexpr size_t N = Shape::N + 1;
    /// Point type of the base shape.
    typedef typename PointType<T,Shape::N>::point_t base_point_t;
    
    /// Cross-section of this frustum at the `h=1` plane.
    Shape base;
    /// Height range spanned by this frustum.
    Rect<T,1> height;
    
    
    /**
     * @brief Construct a pyramid with its tip at the origin
     * and a default-constructed cross-sectional base at h=1.
     */
    Frustum():height(0,1) {}
    
    /**
     * @brief Construct a new Frustum having `base` as a cross section,
     * and spanning heights between `h0` and `h1`.
     */
    Frustum(const Shape& base, T h0, T h1):
                base(base), 
                height(std::min(h0,h1), std::max(h0,h1)) {}
    
    /**
     * @brief Construct a new Frustum having `base` as a cross section,
     * and spanning the height range `h`.
     */
    Frustum(const Shape& base, const Rect<T,1>& h):
        base(base),
        height(h) {}
    
    static constexpr bool admits_cusps() { return true; }
    
    bool operator==(const Frustum& other) const {
        return base == other.base && height == other.height;
    }
    
    /**
     * Frustum-point intersection test.
     * 
     * @param p A point.
     * @return `true` if `p` is on or inside this frustum; `false` otherwise.
     */
    bool contains(Vec<T,N> p) const requires RegionObject<Shape> {
        if (not clipped_height().contains(p[N-1])) return false;
        if (p[N-1] == 0) return p == (T)0;
        p /= p[N-1];
        return base.contains(p.template resized<N-1>());
    }
    
    /// Ray-shape intersection.
    Rect<T,1> intersect(const Ray<T,N>& ray) const requires RayIntersectableObject<Shape> {
        const Rect<T,1> h_range = clipped_height();
        const T h_origin = ray.origin[N - 1];
        const T h_direction = ray.direction[N - 1];
        Rect<T,1> slab_range = h_range.intersect(
            Ray<T,1>(h_origin, h_direction));
        if (slab_range.is_empty()) return Rect<T,1>::empty;

        const base_point_t base_origin = ray.origin.template resized<N - 1>();
        const base_point_t base_direction = ray.direction.template resized<N - 1>();

        if (h_direction == 0) {
            if (h_origin != 0) {
                // At constant nonzero height, homogeneous projection is affine,
                // so the base-shape ray parameter is the original ray parameter.
                return base.intersect(Ray<T,N - 1>(
                    base_origin / h_origin,
                    base_direction / h_origin)) & slab_range;
            }

            // The whole ray lies in the h=0 plane.  A (convex) frustum contains
            // only its tip in this plane.
            if (not h_range.contains((T)0)) return Rect<T,1>::empty;
            index_t pivot = 0;
            for (index_t i = 1; i < N - 1; ++i) {
                if (std::abs(base_direction[i]) > std::abs(base_direction[pivot])) pivot = i;
            }
            if (base_direction[pivot] == 0) {
                return base_origin == (T)0 ? Rect<T,1>::full : Rect<T,1>::empty;
            }
            T s = -base_origin[pivot] / base_direction[pivot];
            return ray.at_multiple(s) == (T)0
                ? Rect<T,1>(s, s)
                : Rect<T,1>::empty;
        }

        if (h_range.lo == 0 and h_range.hi == 0) {
            // A zero-height frustum is just its tip.
            const T tip_s = -h_origin / h_direction;
            return ray.at_multiple(tip_s) == (T)0
                ? Rect<T,1>(tip_s, tip_s)
                : Rect<T,1>::empty;
        }

        // Put a reference point on h=1. Under homogeneous projection to that
        // plane, the original straight ray remains a straight line:
        //
        //   x(s) / h(s) = x0 + t(s) v0
        //   t(s) = (s - s0) / h(s),  h(s0) = 1.
        //
        // This projective parameterization avoids the several nearly-dependent
        // orthogonal projections used by the old implementation.
        const T s0 = ((T)1 - h_origin) / h_direction;
        const base_point_t x0 = base_origin + s0 * base_direction;
        const base_point_t v0 = base_direction - h_direction * x0;

        if (v0 == (T)0) {
            // The ray passes through the tip: its homogeneous projection is a point.
            return base.contains(x0) ? slab_range : Rect<T,1>::empty;
        }

        Rect<T,1> t_range = base.intersect(Ray<T,N - 1>(x0, v0));
        if (t_range.is_empty()) return Rect<T,1>::empty;

        // Restrict the projective parameter before mapping back. The map has a
        // pole at h=0, but clipped_height() guarantees that the retained slab is
        // on only one side of that pole. Its zero-height endpoint maps to the
        // appropriate infinity.
        const T h_at_s_lo = h_direction > 0 ? h_range.lo : h_range.hi;
        const T h_at_s_hi = h_direction > 0 ? h_range.hi : h_range.lo;
        const T neg_inf = -std::numeric_limits<T>::infinity();
        const T pos_inf =  std::numeric_limits<T>::infinity();
        const T t_lo = h_at_s_lo == 0
            ? neg_inf
            : (h_at_s_lo - (T)1) / (h_direction * h_at_s_lo);
        const T t_hi = h_at_s_hi == 0
            ? pos_inf
            : (h_at_s_hi - (T)1) / (h_direction * h_at_s_hi);
        t_range &= Rect<T,1>(t_lo, t_hi);
        if (t_range.is_empty()) return Rect<T,1>::empty;

        const T tip_s = -h_origin / h_direction;
        auto ray_parameter = [&](T t) {
            return std::isinf(t) ? tip_s : s0 + t / ((T)1 - h_direction * t);
        };
        return Rect<T,1>(ray_parameter(t_range.lo), ray_parameter(t_range.hi));
    }

    Vec<T,N> convex_support(Vec<T,N> d) const requires ConvexObject<Shape> {
        Rect<T,1> h = clipped_height();
        T sign = h.lo < 0 ? -1 : 1; // ← the shape is flipped below the origin
        
        // find the appropriate direction to check in shape-space.
        base_point_t d_s = sign * d.template resized<N-1>();
        if (d_s == (T)0) d_s[0] = 1;  // ← pick an arbitrary direction
        
        // find the relevant point on the base shape, and place it on the h=1 plane
        Vec<T,N> p(base.convex_support(d_s), 1);
        
        // top or bottom face?
        // the support point definitely lies on the line passing through 
        // the origin and `p`. which of the two vertices on that line is the 
        // point? the one in the +p direction, or the -p direction?
        T z = (p.dot(d) > 0) ? h.hi : h.lo;
        
        // rescale the point appropriately:
        return z * p;
    }
    
    Rect<T,N> bounds() const requires BoundedObject<Shape> {
        Rect<T,N-1> b0 = base.bounds();
        Rect<T,1> h    = clipped_height();
        // make two rects for the top and bottom faces;
        // union them; extrude by the height:
        return ((b0 * h.lo) | (b0 * h.hi)) * h; // purdy!!
    }
    
    template <ConvexObject S>
    requires ConvexObject<Shape> and (Shape::N == N) and std::same_as<typename S::elem_t, T>
    bool intersects(const S& other) const {
        return geom::intersects(
            as_any_convex(*this),
            as_any_convex(other)
        );
    }
    
    /// Return the height range of this Frustum after it has been clipped by the origin.
    inline Rect<T,1> clipped_height() const {
        const Rect<T,1>& h = height; // shorthand
        return Rect<T,1>((h.lo < 0 and h.hi > 0) ? 0 : h.lo, h.hi);
    }

}; // class Frustum

/** @addtogroup shape
 *  @{
 */

/**
 * @brief Convenience function to raise the shape `s` into a frustum between 
 * heights `h0` and `h1`, by wrapping `s` in the `Frustum` template.
 *
 * @related Frustum
 */
template <typename Shape>
inline Frustum<Shape> frustum(
        const Shape& s,
        typename Shape::elem_t h0,
        typename Shape::elem_t h1)
{
    return Frustum<Shape>(s, h0, h1);
}

/**
 * @brief Convenience typedef for oriented, rectangular, N-dimensional Frustums
 * @related Oriented
 * @related Frustum
 */
template <typename T, index_t N>
using ViewFrustum = Transformed< Frustum< Rect<T,N-1> > >;

/**
 * @brief Convenience typedef for a 3D cone with its tip at the origin.
 * 
 * Formed from a frustum with a circular base.
 * 
 * @related Frustum
 * @tparam T Coordinate type
 */
template <typename T>
using Cone = Frustum<Circle<T>>;

/**
 * @brief Convenience typedef for an oriented cone.
 * 
 * @related Cone
 * @related Transformed
 * 
 * @tparam T Coordinate type
 */
template <typename T>
using TransformedCone = Transformed<Cone<T>>;

/// @} // group shape


template <typename Shape, typename H>
struct Digest<Frustum<Shape>, H> {
    H operator()(const Frustum<Shape> &s) const {
        H nonce = geom::truncated_constant<H>(0x28211b7d8ba5f09b, 0xbcd99a70bc779985);
        return geom::hash_many<H>(nonce, s.base, s.height);
    }
};

#ifdef GEOMC_USE_STREAMS

template <typename Shape>
std::ostream& operator<<(std::ostream& os, const Frustum<Shape>& f) {
    os << "Frustum(" << f.base << ", " << f.height << ")";
    return os;
}

#endif

} // namespace geom


template <typename Shape>
struct std::hash<geom::Frustum<Shape>> {
    size_t operator()(const geom::Frustum<Shape> &s) const {
        return geom::hash<geom::Frustum<Shape>, size_t>(s);
    }
};
