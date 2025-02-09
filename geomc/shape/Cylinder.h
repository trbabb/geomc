#pragma once
/* 
 * File:   Cylinder.h
 * Author: tbabb
 *
 * Created on May 19, 2014, 10:58 PM
 */
#include <geomc/linalg/Similarity.h>
#include <geomc/shape/ShapeTypes.h>
#include <geomc/shape/Shape.h>
#include <geomc/function/Utils.h>

namespace geom {

/** 
 * @ingroup shape
 * @brief An N-dimensional cylinder, given by its radius and endpoints.
 * 
 * Represents an extrusion of an N-2 sphere. In other words, an extrusion of 
 * a disk in 3D; an extrusion of a sphere in 4D; and an extrusion of a line 
 * segment in 2D.
 * 
 * If `p0 == p1`, the behavior is undefined.
 */
template <typename T, index_t N>
class Cylinder: 
    public Convex          <T,N,Cylinder<T,N>>,
    public RayIntersectable<T,N,Cylinder<T,N>>,
    public Projectable     <T,N,Cylinder<T,N>>
{
    public:
    /// Axis endpoint.
    Vec<T,N> p0;
    /// Axis endpoint.
    Vec<T,N> p1;
    /// Cylinder radius.
    T radius;
    
    /// Construct a cylinder of radius and length 1, with axis along X+.
    Cylinder(): radius(1) { p1[0] = 1; }
    
    /// Construct a cylinder of radius `r` and length 1, with axis along X+.
    Cylinder(T r): radius(r) { p1[0] = 1; }
    
    /**
     * @brief Construct a cylinder with arbitrary radius and endpoints.
     * @param p0 An endpoint of the cylinder axis.
     * @param p1 An endpoint of the cylinder axis.
     * @param radius Radius of cylinder.
     */
    Cylinder(const Vec<T,N>& p0, const Vec<T,N>& p1, T radius):
        p0(p0),
        p1(p1),
        radius(radius) {}
    
    static constexpr bool admits_cusps() { return false; }
    
    /// Shape-point intersection test.
    bool contains(Vec<T,N> p) const {
        Vec<T,N> a = p1 - p0;
        Vec<T,N> b = p  - p0;
        T b2 = b.mag2();
        T a2 = a.mag2();
        T d  = a.dot(b);
        T x2 = d * d / a2;
        bool in_r = b2 - x2 <= radius * radius;
        bool in_caps = x2 < a2 and d >= 0;
        return in_r and in_caps;
    }
    
    bool operator==(const Cylinder<T,N>& other) const {
        return p0 == other.p0 and p1 == other.p1 and radius == other.radius;
    }
    
    /// Signed distance function.
    T sdf(Vec<T,N> p) const {
        // approach: put p into the r, a_hat basis,
        // where the nearest cap is at the origin.
        // axis of the cylinder:
        Vec<T,N> a = p1 - p0;
        // vec from "base" of cylinder to p
        Vec<T,N> b = p - p0;
        T b2 = b.mag2(); // squared length of b (hypoteneuse)
        T a2 = a.mag2(); // squared length of axis (side 1)
        T d  = a.dot(b);
        // fractional distance along axis from 0 to 1
        T s  = d / a2;
        // square of distance from the base to the projected axis point
        T x2 = d * s;
        // signed distance of p to the surface of the cylinder
        T r = std::sqrt(b2 - x2) - radius;
        // signed fractional distance along the axis to the nearest cap
        s = std::max(-s, (s - 1));
        // sdf is negative (inside shape) iff (s, r) in the negative quadrant
        T sign = (s < 0 and r < 0) ? -1 : 1;
        // clamped coordinates, to project orthogonally to the wall/cap
        s = std::max(s, (T)0);
        r = std::max(r, (T)0);
        // squared axis-parallel distance to cap
        x2 = s * s * a2;
        // signed distance to the surface point
        return sign * std::sqrt(r * r + x2);
    }
    
    Vec<T,N> project(Vec<T,N> p) const {
        // axis of the cylinder
        Vec<T,N> a = p1 - p0;
        // vec from "base" of cylinder to p
        Vec<T,N> b = p - p0;
        T b2 = b.mag2(); // squared length of b (hypoteneuse)
        T a2 = a.mag2(); // squared length of axis (side 1)
        T d  = a.dot(b);
        // fractional distance along axis from 0 to 1
        T s  = d / a2;
        // square of distance from the base to the projected axis point
        T x2 = d * s;
        // distance of p to the axis
        T r_dist = std::sqrt(b2 - x2);
        // unit length radial vector
        Vec<T,N> r = (p - p0 - s * a) / r_dist;
        // adjust `r` and `s` to lie on the cylinder boundary
        if (s < 0 or s > 1) {
            // p projects onto the cap
            r_dist = std::min(r_dist, radius);
            s = s < 0 ? 0 : 1;
        } else if (r_dist < radius and 
                std::sqrt(a2) * std::min(s, 1 - s) < radius - r_dist) {
            // p is interior and the cap is closer than the wall
            s = s < 0.5 ? 0 : 1;
        } else {
            // p projects onto the cylinder wall
            r_dist = radius;
        }
        return p0 + (s * a) + (r * r_dist);
    }
    
    Vec<T,N> normal(Vec<T,N> p) const {
        // axis of the cylinder
        auto a = p1 - p0;
        // vec from "base" of cylinder to p
        auto b = p - p0;
        T b2 = b.mag2(); // squared length of b (hypoteneuse)
        T a2 = a.mag2(); // squared length of axis (side 1)
        T d  = a.dot(b);
        // fractional distance along axis from 0 to 1
        T s  = d / a2;
        // square of distance from the base to the projected axis point
        T x2 = d * s;
        // distance of p to the axis
        T r_dist = std::sqrt(b2 - x2);
        // unit length radial vector
        auto r = (p - p0 - s * a) / r_dist;
        // choose the outward direction
        Vec<T,N> n;
        if (s < 0 or s > 1) {
            // p projects onto the cap
            if (r_dist > radius) {
                // p projects onto the rim of the cap
                n = p - (r * radius + p0);
            } else {
                // p projects to the cap face; normal is axial
                n = s < 0 ? -a : a;
            }
        } else if (r_dist < radius and std::sqrt(a2) * std::min(s, 1 - s) < radius - r_dist) {
            // p is interior and the cap is closer than the wall
            n = s < 0.5 ? -a : a;
        } else {
            // p projects onto the cylinder wall
            n = r;
        }
        return n.unit();
    }
    
    /**
     * @return An axis-aligned bounding box completely containing this
     * cylinder.
     */
    Rect<T,N> bounds() const {
        // construct two bounding boxes, one for each disk cap.
        // we do this by finding the extent of the cap along each axis.
        Vec<T,N> a  = p1 - p0;
        // square of unit normal's components
        // aka square of projection of normal onto axis of interest
        Vec<T,N> n2 = (a * a) / a.mag2(); 
        // solve for the remaining side of the triangle
        // which happens to be the projection of a perpendicular vector
        Vec<T,N> c0 = radius * std::sqrt((Vec<T,N>::ones - n2));
        // construct the disk bounds at the end cap positions
        Rect<T,N> b0(p0 - c0, p0 + c0);
        Rect<T,N> b1(p1 - c0, p1 + c0);
        // box union
        return b0 | b1;
    }
    
    Vec<T,N> convex_support(Vec<T,N> d) const {
        Vec<T,N> a = p1 - p0;
        Vec<T,N> perp = (d - d.project_on(a)).unit() * radius;
        return (d.dot(a) > 0 ? p1 : p0) + perp;
    }
    
    /// Measure the interior (volume) of the cylinder.
    T measure_interior() const {
        return measure_ball_interior(N - 1, radius) * (p1 - p0).mag();
    }
    
    /// Measure the boundary (surface area) of the cylinder.
    T measure_boundary() const {
        return (
            // shaft area
            measure_sphere_boundary(N - 1, radius) * (p1 - p0).mag() +
            // two endcaps
            2 * measure_ball_interior(N - 1, radius)
        );
    }
    
    /// Ray/shape intersection test.
    Rect<T,1> intersect(const Ray<T,N>& ray) const {
        Vec<T,N> a = p1 - p0;
        Vec<T,N> b = ray.origin - p0;
        T a2 = a.mag2();
        T b2 = b.mag2();
        T d  = a.dot(b);
        T s  = d / a2;  // fractional projection of ray origin onto axis
        T x2 = d * s;
        T m  = ray.direction.dot(a);
        
        // construct quadratic coeffs
        T k0 = ray.direction.mag2() - m * m / a2;
        T k1 = 2 * (b.dot(ray.direction) - d * m / a2);
        T k2 = b2 - x2 - radius * radius;
        
        // setup for cylinder intersection
        T roots[2];
        Rect<T,1> interval;
        
        // construct an interval holding the ray's overlap with the infinite cylinder
        if (k0 != 0 and quadratic_solve(roots, k0, k1, k2)) {
            interval = Rect<T,1>::spanning_corners(roots[0], roots[1]);
        } else if (k2 < 0) {
            // if there is no intersection, then it might be the case that 
            // the ray points exactly down the interior of the cylinder, which
            // can only happen if the ray origin is inside too (k2 < 0). also
            // due to precision this manifests sometimes when there is one root:
            // when v and a are aligned then k0 is 0, and the equation can take the form
            // `0x^2 + (small)x + c = 0`; this should be treated the same as 'no solution'.
            // (Tangent rays with a single root will have k0 ≠ 0 and should be intersected
            // normally).
            interval = Rect<T,1>::full;
        } else {
            // no intersection
            return Rect<T,1>::empty;
        }
        
        // trace the endcaps
        if (m != 0) {
            // trace the two cap planes
            T z0 = d; // a.dot(ray.origin - p0);
            T z1 =       a.dot(ray.origin - p1);
            T s0 = -z0 / m;
            T s1 = -z1 / m;
            // intersect the slab with the cylinder
            interval &= Rect<T,1>::spanning_corners(s0, s1);
        } else if (s < 0 or s > 1) {
            // ray is parallel to the slab and outside of it
            return Rect<T,1>::empty;
        }
        
        return interval;
    }
    
};

/// @addtogroup shape
/// @{

/// @brief Transform a cylinder by a similarity transform.
/// @related Cylinder
/// @related Similarity
template <typename T, index_t N>
Cylinder<T,N> operator*(const Similarity<T,N>& xf, const Cylinder<T,N>& c) {
    return Cylinder<T,N>(
        xf * c.p0,
        xf * c.p1,
        c.radius * xf.sx
    );
}

/// @brief Inverse-transform a cylinder by a similarity transform.
/// @related Cylinder
/// @related Similarity
template <typename T, index_t N>
Cylinder<T,N> operator/(const Cylinder<T,N>& c, const Similarity<T,N>& xf) {
    return Cylinder<T,N>(
        xf / c.p0,
        xf / c.p1,
        c.radius / xf.sx
    );
}

/// @brief Transform a cylinder by an isometry.
/// @related Cylinder
/// @related Isometry
template <typename T, index_t N>
Cylinder<T,N> operator*(const Isometry<T,N>& xf, const Cylinder<T,N>& c) {
    return Cylinder<T,N>(
        xf * c.p0,
        xf * c.p1,
        c.radius
    );
}

/// @brief Inverse-transform a cylinder by an isometry.
/// @related Cylinder
/// @related Isometry
template <typename T, index_t N>
Cylinder<T,N> operator/(const Cylinder<T,N>& c, const Isometry<T,N>& xf) {
    return Cylinder<T,N>(
        xf / c.p0,
        xf / c.p1,
        c.radius
    );
}

/// @} // addtogroup shape

template <typename T, index_t N, typename H>
struct Digest<Cylinder<T,N>,H> {
    H operator()(const Cylinder<T,N> &c) const {
        H nonce = geom::truncated_constant<H>(0x77f68ad97f8281e6, 0xb0ffe271f8a4f531);
        return geom::hash_many<H>(nonce, c.p0, c.p1, c.radius);
    }
};

    
} // namespace geom


template <typename T, index_t N>
struct std::hash<geom::Cylinder<T,N>> {
    size_t operator()(const geom::Cylinder<T,N> &v) const {
        return geom::hash<geom::Cylinder<T,N>, size_t>(v);
    }
};
