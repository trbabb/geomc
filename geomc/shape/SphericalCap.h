#pragma once

#include <numbers>
#include <geomc/linalg/Vec.h>
#include <geomc/shape/Shape.h>
#include <geomc/shape/Sphere.h>

namespace geom {

/**
 * @brief The surface of a sphere within a certain angle from its pole.
 * 
 * The cap is defined by a half-angle, which is the angle between the pole and
 * the edge of the cap. The pole is the final axis in N-dimensional space.
 *
 * The cap is a surface, and has no interior. It is concave towards the origin.
 * 
 * In 2D, the cap is an arc centered on the Y+ axis. In 3D, it's a cap centered on Z+.
 * 
 * It is valid for the cap to extend beyond the equator. The half angle is
 * canonically between 0 and π.
 * 
 * To orient or scale the cap, wrap it with a Similar<Shape>.
 * 
 * @ingroup shape
 */
template <typename T, index_t N>
class SphericalCap : public Dimensional<T,N> {
public:
    using typename Dimensional<T,N>::point_t;
    
    static_assert(N > 1, "SphericalCap must be at least 2D.");
    
    /// Angle between the polar axis and the edge of the cap, in radians.
    T half_angle_radians;
    
    /// Construct a spherical cap with the given half-angle.
    SphericalCap(T radians):
        half_angle_radians(radians) {}
        
    static constexpr bool admits_cusps() { return true; }
    
    /// Shape equality.
    bool operator==(const SphericalCap& other) const {
        return half_angle_radians == other.half_angle_radians;
    }
    
    /// Axis aligned bounding box.
    Rect<T,N> bounds() const {
        T s = std::sin(half_angle_radians);
        T c = std::cos(half_angle_radians);
        return Rect<T,N-1>(-s, s) * Rect<T,1>(c, 1);
    }
    
    /// Signed distance function.
    T sdf(point_t p) const {
        if (p.is_zero()) return 1;
        if (is_inside_angle(p)) {
            // distance to the sphere surface
            return p.mag() - 1;
        } else {
            // distance to the boundary of the cap
            VecType<T,N-1> p_proj = p.template resized<N-1>();
            T s = std::sin(half_angle_radians);
            T c = std::cos(half_angle_radians);
            T z = p[N - 1];
            T x = mag(p_proj);
            T dx = x - s;
            T dz = z - c;
            return std::sqrt(dx * dx + dz * dz);
        }
    }
    
    /**
     * @brief Point containment test. Always false, because the shape is a surface.
     *
     * This is implemented to conform to the SdfObject concept, since a signed distance
     * function trivially implies the ability to test for point containment. The signed
     * distance function never returns a value less than zero, so this method always
     * returns false.
     */
    constexpr bool contains(point_t p) const {
        return false;
    }
    
    /// Project the point `p` orthogonally onto the surface of the cap.
    point_t project(point_t p) const {
        if (p.is_zero()) {
            // project the origin to the pole
            return point_t(VecType<T,N-1>((T)0), 1);
        } else if (is_inside_angle(p)) {
            // project to sphere surface
            return p.unit();
        } else {
            // project to the boundary of the cap
            VecType<T,N-1> p_proj = p.template resized<N-1>();
            T s = std::sin(half_angle_radians);
            T c = std::cos(half_angle_radians);
            if constexpr (N > 2) {
                return point_t(p_proj.with_length(s), c);
            } else {
                return point_t(s, c);
            }
        }
    }
    
    /// Outward-facing direction.
    point_t normal(point_t p) const {
        if (p.is_zero()) {
            // normal to the pole
            return point_t(VecType<T,N-1>((T)0), -1);
        } else if (is_inside_angle(p)) {
            // normal to the sphere surface
            T sign = p.mag2() > 1 ? 1 : -1;
            return sign * p.unit();
        } else {
            // normal to the boundary of the cap.
            // project p to the boundary:
            VecType<T,N-1> p_proj = p.template resized<N-1>();
            T s = std::sin(half_angle_radians);
            T c = std::cos(half_angle_radians);
            point_t p_boundary = {p_proj.with_length(s), c};
            return (p - p_boundary).unit();
        }
    }
    
    /**
     * @brief Clamp the coordinates of `p` to lie within this shape.
     *
     * This is the same as project(), because the shape is a surface and has no interior.
     */
    point_t clip(point_t p) const {
        return project(p);
    }
    
    /**
     * @brief Intersect the cap with a ray.
     * 
     * The ray may intersect the cap at 0, 1, or 2 points. In the case of 2 points,
     * the ray parameters will be the two endpoints of the interval. In the case of 1 point,
     * the endpoints of the interval will be identical. In the case of 0 points, the interval is
     * empty.
     *
     * Note that, unlike with other shapes, the interval between the two endpoints is not
     * part of the cap.
     */
    Rect<T,1> intersect(const Ray<T,N>& ray) const {
        Rect<T,1> sphere_interval = Sphere<T,N>().intersect(ray);
        if (sphere_interval.is_empty()) return sphere_interval;
        T ab[2] = {sphere_interval.lo, sphere_interval.hi};
        Rect<T,1> out_interval;
        for (index_t i = 0; i < 2; ++i) {
            Vec<T,N> p = ray * ab[i];
            if (is_inside_angle(p)) {
                out_interval |= ab[i];
            }
        }
        return out_interval;
    }
    
    /**
     * @brief Measure the boundary (surface area) of the shape.
     *
     * For a 2D spherical cap, this is the arc length; the same as its full angle.
     * For 3D, this is the area of the cap; the same as its solid angle.
     *
     * Higher dimensions are not yet available, but may be provided in the future.
     */
    T measure_boundary() const requires (N == 2 or N == 3) {
        // todo: measure higher dimensions (requires some ugly math)
        if constexpr (N == 2) {
            return 2 * half_angle_radians;
        } else if constexpr (N == 3) {
            return 2 * std::numbers::pi_v<T> * (1 - std::cos(half_angle_radians));
        }
    }
    
    /**
     * @brief Test whether the point `p` is inside the cap's angle.
     *
     * Returns `true` if `p` orthogonally projects to the cap's surface;
     * `false` if it projects to the cap's edge. 
     */
    bool is_inside_angle(point_t p) const {
        T m2 = p.mag2();
        T z  = p[N - 1];
        if (m2 == 0) return true;
        T t = geom::positive_mod<T>(half_angle_radians / std::numbers::pi_v<T>, 2);
        // half-turns from the axis, in [-0.5, 0.5]
        // 0 is the equatorial plane, -0.5 is +pole, 0.5 is -pole
        T u = 0.5 - std::abs(1 - t);
        // split into two cases for better numerical precision
        if (std::abs(u) < 0.25){
            // point is closer to the equatorial plane;
            // test based on the height above the equatorial plane,
            // where angle mostly affects height
            T c = std::cos(half_angle_radians);
            T s_cap = c < 0 ? -1 : 1;
            T s_pt  = z < 0 ? -1 : 1;
            return s_pt * z * z / m2 >= s_cap * c * c;
        } else {
            // point is closer to the axis;
            // test based on the distance from the axis,
            // where angle mostly affects radius
            T s  = std::sin(half_angle_radians);
            VecType<T,N-1> p_base = p.template resized<N-1>();
            T k = mag2(p_base) / m2;
            if (u < 0) {
                // angle is above the equator; area nearest the axis is included.
                // everything below the equator is excluded.
                return k <= s * s and z > 0;
            } else {
                // angle is below the equator; area nearest the axis is excluded.
                // everything above the equator is included.
                return k >= s * s or z > 0;
            }
        }
    }
    
};

} // namespace geom
