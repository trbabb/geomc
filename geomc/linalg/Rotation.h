#pragma once

#include <cmath>
#include <numbers>
#include <geomc/linalg/Quaternion.h>
#include <geomc/linalg/AffineTransform.h>

namespace geom {

// nb: we use multiplication to denote rotation composition, because
//   it is in general not commutative, and because it is a group operation.
//   it would be misleading to use '+' for a non-commutative operation.

/**
 * @ingroup linalg
 * @brief A rotation in N-dimensional space.
 * 
 * Currently 2D and 3D are implemented. See Rotation<T,2> and Rotation<T,3>.
 *
 * Rotations meet the Transform concept.
 *
 * For transforms which include a translation, see Isometry.
 *
 * For transforms which include a translation and a scaling, see Similarity.
 * 
 * For nonuniform scaling or skew transforms, see AffineTransform.
 *
 * Rotations are composed like transforms, with multiplication on the left:
 *
 *     Rotation<T,N> r1, r2;
 *     Rotation<T,N> r3 = r2 * r1; // rotation which applies r2, then r1
 *     Vec<T,N> v = r3 * r1 * v0;  // apply r1 to v0, then r3 to the result
 * 
 * Rotations can be inverted with the `/` operator:
 *
 *     Rotation<T,N> r1, r2;
 *     Rotation<T,N> r3 = r2 / r1; // rotation which takes r1 to r2
 *     Vec<T,N> v = v0 / r3;       // apply the inverse of r3 to v0
 *
 */
template <typename T, index_t N>
class Rotation {}; // todo: implement N-blades

/**
 * @ingroup linalg
 * @brief 2D rotation.
 *
 * See `Rotation` for a general description of rotations.
 */
template <typename T>
class Rotation<T,2> : public Dimensional<T,2> {
public:
    
    T radians;
    
    Rotation():radians(0) {}
    Rotation(T radians):radians(radians) {}
    
    static Rotation<T,2> align_vectors(const Vec<T,2>& dir, const Vec<T,2>& align_with) {
        // numerically stable formula from https://scicomp.stackexchange.com/a/27769
        Vec<T,2> a = align_with.unit();
        Vec<T,2> b = dir.unit();
        T d_a = 2. * std::atan2((a + b).mag(), (a - b).mag());
        return Rotation(d_a);
    }
    
    /// Represent as an affine transform.
    AffineTransform<T,2> transform() const {
        return geom::rotation(radians);
    }
    
    /// Cast to affine transform.
    operator AffineTransform<T,2>() const {
        return transform();
    }
    
    /// Cast to complex number.
    operator std::complex<T>() const {
        return std::polar((T)1, radians);
    }
    
    /// Cast the underlying coordinate type.
    template <typename U>
    explicit operator Rotation<U,2>() const {
        return Rotation<U,2>(static_cast<U>(radians));
    }
    
    /// Cast to a 3D rotation.
    operator Rotation<T,3>() const {
        return Rotation<T,3>({0,0,1}, radians);
    }
    
    /// Cast to a 3D rotation and change the coordinate type.
    template <typename U>
    explicit operator Rotation<U,3>() const {
        return Rotation<U,3>({0,0,1}, static_cast<U>(radians));
    }
    
    /// Apply to a 2D vector.
    Vec<T,2> operator*(const Vec<T,2>& v) const {
        return v.rotated(radians);
    }
    
    /// Apply to a direction vector (conorming to Transform concept)
    Vec<T,2> apply_direction(const Vec<T,2>& v) const {
        return (*this) * v;
    }
    
    /// Compose rotation.
    Rotation<T,2> operator*(const Rotation<T,2>& other) const {
        return Rotation<T,2>(radians + other.radians);
    }
    
    /// In-place compose rotation.
    Rotation<T,2>& operator*=(const Rotation<T,2>& other) {
        radians += other.radians;
        return *this;
    }
    
    /// Find the rotation that takes `other` to `this`.
    Rotation<T,2> operator/(const Rotation<T,2>& other) const {
        return Rotation<T,2>(radians - other.radians);
    }
    
    /// In-place apply inverse
    Rotation<T,2>& operator/=(const Rotation<T,2>& other) {
        radians -= other.radians;
        return *this;
    }
    
    /// Apply inverse rotation to a direction vector (conforming to Transform concept)
    Vec<T,2> apply_inverse_direction(const Vec<T,2>& v) const {
        return v.rotated(-radians);
    }
    
    /// In-place scaling of a rotation.
    Rotation<T,2>& operator*=(T s) {
        radians *= s;
        return *this;
    }
    
    Rotation<T,2> exp() const {
        return *this;
    }
    
    Rotation<T,2> log() const {
        return *this;
    }
    
    /// Return the angle component of the rotation.
    T angle() const {
        return radians;
    }
    
    /// Return a rotation with the same orientation, but normalized to the rotation
    /// angle range [0, 2π).
    Rotation<T,2> canonical() const {
        return Rotation<T,2>(
            geom::positive_mod<T>(radians, std::numbers::pi_v<T> * 2)
        );
    }
    
    /// Compute the inverse rotation.
    Rotation<T,2> inv() const {
        return Rotation<T,2>(-radians);
    }
};

/** @addtogroup linalg
 *  @{
 */

/// @brief Apply the inverse of a rotation to a vector.
/// @related Rotation
template <typename T>
Vec<T,2> operator/(const Vec<T,2>& v, const Rotation<T,2>& r) {
    return v.rotated(-r.radians);
}

/// @brief Apply a rotation to a ray.
/// @related Rotation
/// @related Ray
template <typename T, index_t N>
Ray<T,N> operator*(const Rotation<T,N>& rot, Ray<T,N> ray) {
    return {rot * ray.direction, rot * ray.origin};
}

/// @brief Apply the inverse of a rotation to a ray.
/// @related Rotation
/// @related Ray
template <typename T, index_t N>
Ray<T,N> operator/(Ray<T,N> ray, const Rotation<T,N>& rot) {
    return {rot / ray.direction, rot / ray.origin};
}

/// @brief Extend a rotation.
/// @related Rotation
template <typename T>
Rotation<T,2> operator*(T s, const Rotation<T,2>& o) {
    return {o.radians * s};
}

/// @brief Extend a rotation.
/// @related Rotation
template <typename T>
Rotation<T,2> operator*(const Rotation<T,2>& o, T s) {
    return {o.radians * s};
}

/// @brief Minimally interpolate two rotations.
/// @related Rotation
template <typename T>
Rotation<T,2> mix(T s, const Rotation<T,2>& a, const Rotation<T,2>& b) {
    T d_a = geom::angle_to(a.radians, b.radians);
    return Rotation<T,2>(a.radians + s * d_a);
}

/**
 * @brief 3D rotation.
 *
 * See `Rotation` for a general description of rotations.
 */
template <typename T>
class Rotation<T,3> : public Dimensional<T,3> {
public:
    
    Quat<T> q;
    
    Rotation():q(0,0,0,1) {}
    Rotation(const Quat<T>& q):q(q) {}
    Rotation(Vec<T,3> axis, T angle):
        q(Quat<T>::rotation_from_axis_angle(axis, angle)) {}
    
    static Rotation<T,3> align_vectors(const Vec<T,3>& dir, const Vec<T,3>& align_with) {
        return {
            Quat<T>::rotation_direction_align(dir, align_with)
        };
    }
    
    /// Represent as an affine transform.
    AffineTransform<T,3> transform() const {
        return geom::rotation(q);
    }
    
    /// Cast to affine transform.
    operator AffineTransform<T,3>() const {
        return transform();
    }
    
    /// Cast to quaternion.
    operator Quat<T>() const {
        return q;
    }
    
    /// Apply to a 3D vector.
    Vec<T,3> operator*(const Vec<T,3>& v) const {
        return q * v;
    }
    
    /// Apply to a direction vector (conforming to Transform concept)
    Vec<T,3> apply_direction(const Vec<T,3>& v) const {
        return (*this) * v;
    }
    
    /// Compose rotation.
    Rotation<T,3> operator*(const Rotation<T,3>& other) const {
        return Rotation<T,3>(q * other.q);
    }
    
    /// In-place compose rotation.
    Rotation<T,3>& operator*=(const Rotation<T,3>& other) {
        q = q * other.q;
        return *this;
    }
    
    /// Find the rotation that takes `other` to `this`.
    Rotation<T,3> operator/(const Rotation<T,3>& other) const {
        return Rotation<T,3>(q * other.q.conj());
    }
    
    /// In-place apply inverse
    Rotation<T,3>& operator/=(const Rotation<T,3>& other) {
        q = q * other.q.conj();
        return *this;
    }
    
    /// Apply inverse rotation to a direction vector (conforming to Transform concept)
    Vec<T,3> apply_inverse_direction(const Vec<T,3>& v) const {
        return q.conj() * v;
    }
    
    /// In-place scaling of a rotation.
    Rotation<T,3>& operator*=(T s) {
        q = q.rotation_scale(s);
        return *this;
    }
    
    Rotation<T,3> exp() const {
        return { q.exp() };
    }
    
    Rotation<T,3> log() const {
        return { q.log() };
    }
    
    /// Return the angle component of the rotation.
    T angle() const {
        return q.angle();
    }
    
    /// Return a rotation with the same orientation, but with axis chosen so the rotation
    /// angle is in the range [0, π].
    Rotation<T,3> canonical() const {
        T sign = q.w < 0 ? -1 : 1;
        return Rotation<T,3>(sign * q.unit());
    }
    
    /// Compute the inverse rotation.
    Rotation<T,3> inv() const {
        return Rotation<T,3>(q.conj());
    }
};

/// @brief Apply the inverse of a rotation to a vector.
/// @related Rotation
template <typename T>
Vec<T,3> operator/(const Vec<T,3>& v, const Rotation<T,3>& r) {
    return r.q.conj() * v;
}

/// @brief Extend a rotation.
/// @related Rotation
template <typename T>
Rotation<T,3> operator*(T s, const Rotation<T,3>& o) {
    return {o.q.rotation_scale(s)};
}

/// @brief Extend a rotation.
/// @related Rotation
template <typename T>
Rotation<T,3> operator*(const Rotation<T,3>& o, T s) {
    return {o.q.rotation_scale(s)};
}

/// @brief Minimally interpolate two rotations.
/// @related Rotation
template <typename T>
Rotation<T,3> mix(T s, const Rotation<T,3>& a, const Rotation<T,3>& b) {
    Quat<T> q = a.q.slerp(b.q, s);
    return Rotation<T,3>(q);
}

/// @} // addtogroup linalg

template <typename T, typename H>
struct Digest<Rotation<T,2>, H> {
    H operator()(const Rotation<T,2>& r) const {
        H nonce = geom::truncated_constant<H>(0x8f4676a597229f40, 0x2f01ddf7f5cb29b9);
        return geom::hash_many(nonce, r.radians);
    }
};

template <typename T, typename H>
struct Digest<Rotation<T,3>, H> {
    H operator()(const Rotation<T,3>& r) const {
        H nonce = geom::truncated_constant<H>(0xe0da0dc698074138, 0x60197ba709a81c40);
        return geom::hash_many(nonce, r.q);
    }
};

} // namespace geom
