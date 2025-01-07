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
 * @brief A rotation in N-dimensional space.
 * 
 * Currently 2D and 3D are implemented.
 */
template <typename T, index_t N>
struct Rotation {}; // todo: implement N-blades

/**
 * @brief 2D rotation.
 */
template <typename T>
struct Rotation<T,2> {
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
    
    /// Apply to a 2D vector.
    Vec<T,2> operator*(const Vec<T,2>& v) const {
        return v.rotated(radians);
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
    Rotation<T,2> normalized() const {
        return Rotation<T,2>(
            geom::positive_mod<T>(radians, std::numbers::pi_v<T> * 2)
        );
    }
    
    /// Compute the inverse rotation.
    Rotation<T,2> inv() const {
        return Rotation<T,2>(-radians);
    }
};

template <typename T>
Rotation<T,2> operator*(T s, const Rotation<T,2>& o) {
    return {o.radians * s};
}

template <typename T>
Rotation<T,2> operator*(const Rotation<T,2>& o, T s) {
    return {o.radians * s};
}

/// Minimally interpolate two rotations.
template <typename T>
Rotation<T,2> mix(T s, const Rotation<T,2>& a, const Rotation<T,2>& b) {
    T d_a = geom::angle_to(a.radians, b.radians);
    return Rotation<T,2>(a.radians + s * d_a);
}

/**
 * @brief 3D rotation.
 */
template <typename T>
struct Rotation<T,3> {
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
    
    /// Return a rotation with the same orientation, but normalized to the rotation
    /// angle range [0, 2π).
    Rotation<T,3> normalized() const {
        return Rotation<T,3>(q.unit());
    }
    
    /// Compute the inverse rotation.
    Rotation<T,3> inv() const {
        return Rotation<T,3>(q.conj());
    }
};

template <typename T>
Rotation<T,3> operator*(T s, const Rotation<T,3>& o) {
    return {o.q.rotation_scale(s)};
}

template <typename T>
Rotation<T,3> operator*(const Rotation<T,3>& o, T s) {
    return {o.q.rotation_scale(s)};
}

/// Minimally interpolate two rotations.
template <typename T>
Rotation<T,3> mix(T s, const Rotation<T,3>& a, const Rotation<T,3>& b) {
    Quat<T> q = a.q.slerp(b.q, s);
    return Rotation<T,3>(q);
}

} // namespace geom
