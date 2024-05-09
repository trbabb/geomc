#pragma once

#include <climits>
#include <geomc/shape/Rect.h>
#include <geomc/shape/Sphere.h>

namespace geom {

// todo: point projection is unsolved
//   - can't project point to base plane and recurse
//   - example case: consider a shape far from origin
//     - pt will project to the wall directly "underneath" it, in the limit
//     - in general, projection vector field is "stretched" away from origin by this effect

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
 * templated type alias for `Oriented<Frustum<Rect>>`. 
 */
template <typename Shape>
class Frustum: 
    public Convex          <typename Shape::elem_t, Shape::N+1, Frustum<Shape>>,
    public RayIntersectable<typename Shape::elem_t, Shape::N+1, Frustum<Shape>>
{
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
    
    
    bool operator==(const Frustum& other) const {
        return base == other.base && height == other.height;
    }
    
    /**
     * Frustum-point intersection test.
     * 
     * @param p A point.
     * @return `true` if `p` is on or inside this frustum; `false` otherwise.
     */
    bool contains(Vec<T,N> p) const {
        if (not clipped_height().contains(p[N-1])) return false;
        if (p[N-1] == 0) return p == (T)0;
        p /= p[N-1];
        return base.contains(p.template resized<N-1>());
    }
    
    /// Ray-shape intersection.
    Rect<T,1> intersect(const Ray<T,N>& ray) const {
        // todo: check N=2 case
        // s-values for the top and bottom slabs
        Rect<T,1> slab_range = clipped_height().intersect(
            Ray<T,1>(ray.origin[N-1], ray.direction[N-1]));
        
        // we want to project the ray to the h=1 plane. the family of rays that
        // homogeneously project to a common line on h=1 belong to the plane spanned by 
        // `o` (the ray origin) and `v`, so we will work in that subspace. we want to choose
        // an `o` that is close to the origin for stability, so we ⟂ project Z+ to the (o,v)
        // subspace.
        
        // first make our ray subspace basis orthogonal:
        Vec<T,N> o_ortho = ray.origin - ray.origin.project_on(ray.direction);
        Vec<T,N> z_plus;
        z_plus[N-1] = 1;
        // project the Z+ axis to the ray subspace. (the point on the projected ray
        // which is closest to the origin is a multiple of this point)
        // todo: this can probably be optimized:
        Vec<T,N> o_base = z_plus.project_on(ray.direction) + z_plus.project_on(o_ortho);
        // project V to the h=0 plane
        Vec<T,N> v_base = ray.direction - ray.direction.project_on(o_base);
        
        // handle degenerate cases
        if (o_base[N-1] == 0) {
            // the ray lies in the h=0 plane
            if (o_base == (T)0) [[unlikely]] {
                // the ray passes through the origin, and intersects the
                // frustum tip, if it's included in the height range
                T s = ray.direction.fraction_on(ray.origin);
                return Rect<T,1>(s, s) & slab_range;
            } else {
                // the ray does not pass through the origin
                return Rect<T,1>::empty;
            }
        } else if (v_base == (T)0) {
            // ray passes through the origin. the frustum intersection is full iff the ray is 
            // inside the the frustum; else it is empty.
            return base.contains(o_base.template resized<N-1>())
                ? slab_range         // full range  & slab_range
                : Rect<T,1>::empty;  // empty range & slab_range
        }
        // put `o` onto the h=1 plane
        o_base /= o_base[N-1];
        
        // form the intersection in the base plane
        Ray<T,N-1> ray_base = Ray<T,N-1>(
            o_base.template resized<N-1>(),
            v_base.template resized<N-1>()
        );
        Rect<T,1> s = base.intersect(ray_base);
        
        // short circuit: if the projected ray does not intersect
        // the base shape, then it does not intersect the frustum.
        if (s.is_empty()) {
            return Rect<T,1>::empty;
        }
        
        /* find the hit point `p` on the h=0 plane; the final hit point `q` will lie on some 
           multiple of `p`. The trick is to write `p` and `q` as a points in the (o,v) basis,
           and compute the multiple of `v` that takes us to a multiple of `p`.
           We are solving this geometrical problem:
           
               origin
                  •        ⎫      ⎫
                 ╱ ╲       │      ⎬ p_y
                ╱   • p    ⎬ 1.0  ⎭
               ╱     ·     │
            o •───────•──> ⎭
                      q  v
              
              ╰──────────╯ 1.0
                  ╰───╯    k * p_x
              ╰───╯        o_x
              ╰───────╯    s   <<< want to find this
        */
        Vec<T,N> lo {ray_base.at_multiple(s.lo), 1};
        Vec<T,N> hi {ray_base.at_multiple(s.hi), 1};
        T o_x    = ray.origin.fraction_on(ray.direction);
        T p_lo_y = lo.fraction_on(o_ortho);
        T p_hi_y = hi.fraction_on(o_ortho);
        T p_lo_x = lo.fraction_on(ray.direction);
        T p_hi_x = hi.fraction_on(ray.direction);
        
        // construct the interval of overlap with the infinite frustum.
        s = Rect<T,1>::spanning_corners(
            p_lo_x / p_lo_y - o_x,
            p_hi_x / p_hi_y - o_x
        );
        
        T h0 = ray.origin[N-1] + s.lo * ray.direction[N-1];
        T h1 = ray.origin[N-1] + s.hi * ray.direction[N-1];
        if ((h0 < 0) != (h1 < 0)) {
            // if the two hits are on either side of h=0, then the interval is inverted.
            // (we have already handled the ray through the origin case).
            // which of the two half-infinite intervals should be intersected with the slab?
            s = (slab_range.lo < 0)
                ? Rect<T,1>(std::numeric_limits<T>::lowest(), s.lo)
                : Rect<T,1>(s.hi, std::numeric_limits<T>::max());
        }
        // intersect the slab with the frustum interval.
        return s & slab_range;
    }
    
    Vec<T,N> convex_support(Vec<T,N> d) const {
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
    
    Rect<T,N> bounds() const {
        Rect<T,N-1> b0 = base.bounds();
        Rect<T,1> h    = clipped_height();
        // make two rects for the top and bottom faces;
        // union them; extrude by the height:
        return ((b0 * h.lo) | (b0 * h.hi)) * h; // purdy!!
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
using ViewFrustum = Oriented< Frustum< Rect<T,N-1> > >;

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
 * @related Oriented
 * 
 * @tparam T Coordinate type
 */
template <typename T>
using OrientedCone = Oriented<Cone<T>>;


/** @addtogroup traits
 *  @{
 */

// Frustums inherit concepts
template <typename Shape>
struct implements_shape_concept<Frustum<Shape>, RayIntersectable> : 
    public std::integral_constant<
        bool,
        implements_shape_concept<Shape, RayIntersectable>::value>
{};

template <typename Shape>
struct implements_shape_concept<Frustum<Shape>, Convex> : 
    public std::integral_constant<
        bool,
        implements_shape_concept<Shape, RayIntersectable>::value>
{};


// frustum doesn't yet support sdf or projection.
// template <typename Shape>
// struct implements_shape_concept<Frustum<Shape>, Projectable> : 
//     public std::integral_constant<
//         bool,
//         implements_shape_concept<Shape, Projectable>::value>
// {};

// template <typename Shape>
// struct implements_shape_concept<Frustum<Shape>, SdfEvaluable> : 
//     public std::integral_constant<
//         bool,
//         implements_shape_concept<Shape, Projectable>::value>
// {};

/// @} // addtogroup traits
/// @} // addtogroup shape

template <typename Shape, typename H>
struct Digest<Frustum<Shape>, H> {
    H operator()(const Frustum<Shape> &s) const {
        H nonce = geom::truncated_constant<H>(0x28211b7d8ba5f09b, 0xbcd99a70bc779985);
        return geom::hash_many<H>(nonce, s.base, s.height);
    }
};

} // namespace geom


template <typename Shape>
struct std::hash<geom::Frustum<Shape>> {
    size_t operator()(const geom::Frustum<Shape> &s) const {
        return geom::hash<geom::Frustum<Shape>, size_t>(s);
    }
};
