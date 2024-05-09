/*
 * Quaternion.h
 *
 *  Created on: Mar 30, 2013
 *      Author: tbabb
 */

#ifndef QUATERNION_H_
#define QUATERNION_H_

#include <cmath>
#include <geomc/linalg/Vec.h>

namespace geom {

/** @ingroup linalg 
 *  @brief Quaternion class.
 * 
 * `(x, y, z)` is the vector part, while `w` is the real part. This class differs
 * slightly from convention in that the real part is the last coordinate rather than
 * the first; this scheme was chosen to maintain naming consistency with the
 * rest of the library.
 */
template <typename T>
class Quat : public detail::VecCommon< T, 4, Quat<T> > {
public:
    
    /*******************************
     * Constructors                *
     *******************************/
    
    /// Construct a quaternion with all elements 0
    Quat():detail::VecCommon< T, 4, Quat<T>>() {}
    
    /// Construct a quaternion with elements `(x, y, z, w)`
    Quat(T x, T y, T z, T w) {
            detail::VecBase<T,4>::x = x;
            detail::VecBase<T,4>::y = y;
            detail::VecBase<T,4>::z = z;
            detail::VecBase<T,4>::w = w;
    }
    
    /// Construct a quaternion from the contents of the 4-element array `v`
    Quat(const T v[4]):detail::VecCommon< T, 4, Quat<T>>(v) {}
    
    /// Construct a quaternion with vector part `v` and real part `w`.
    Quat(const Vec<T,3>& v, T w) {
        detail::VecBase<T,4>::x = v.x;
        detail::VecBase<T,4>::y = v.y;
        detail::VecBase<T,4>::z = v.z;
        detail::VecBase<T,4>::w = w;
    }
    
    /// Construct a quaternion from the 4D vector `v`.
    Quat(const Vec<T,4>& v):detail::VecCommon< T, 4, Quat<T>>(v.begin()) {}

    /**
     * Construct a quaternion from a brace-initialization list. (c++11)
     *
     * Example: `Quat<float> q = {0, 0, 0, 1};`
     * @param items A brace-initializer list.
     */
    Quat(const std::initializer_list<T>& items):
            detail::VecCommon< T,4,Quat<T>>(items.begin()) {
        // items.size() not constexpr in c++11  D:<
        static_assert(
            items.size() == 4,
            "Quaternion must be initialized with 4 elements");
    }
    
    /*******************************
     * Static constructors         *
     *******************************/
    
    /**
     * @brief Rotation quaternion from an axis and angle
     * @param axis Axis of rotation (not necessesarily a unit vector)
     * @param angle Angle of rotation in radians
     * @return A unit quaternion representing a rotation of `angle` radians
     * about `axis`.
     */
    static Quat<T> rotation_from_axis_angle(const Vec<T,3> &axis, T angle) {
        if (axis.x == 0 && axis.y == 0 && axis.z == 0) {
            return Quat<T>(axis, 1);
        } else {
            return Quat<T>(axis.unit() * std::sin(angle*0.5), std::cos(angle*0.5));
        }
    }
     
    /**
     * @brief Rotation to align one vector with another.
     * @param v Unit vector to align.
     * @param alignWith Unit direction to align with.
     * @return A quaternion rotating `v` into alignment with `alignWith`.
     */
    static inline Quat<T> rotation_direction_align(
            const Vec<T,3>& v,
            const Vec<T,3>& alignWith)
    {
        Vec<T,3> axis = v ^ alignWith;
        T d = v.dot(alignWith);
        return Quat<T>(axis, d);
    }

    /*******************************
     * Operators                   *
     *******************************/
    
    /// Quaternion multiplication
    friend inline Quat<T> operator*(const Quat<T>& q1, const Quat<T>& q2) {
        return q1.mult(q2);
    }
    
    /**
     * @brief Quaternion conjugation, `q * v * q'`. 
     * 
     * If `q` is a unit quaternion, then `q * v` is a rotation of `v` by `q`.
     */
    friend inline Vec<T,3> operator*(const Quat<T>& q, const Vec<T,3>& v) {
        Quat<T> qv(v,0);
        return (q * qv * q.conj()).imag();
    }
    
    /**
     * @brief Quaternion inverse conjugation, q' * v * q. 
     * 
     * If `q` is a unit quaternion, then `v * q` is a rotation of `v` by `q`<sup>`-1`</sup>.
     */
    friend inline Vec<T,3> operator*(const Vec<T,3>& v, const Quat<T>& q) {
        Quat<T> qv(v,0);
        return (q.conj() * qv * q).imag();
    }
    
    /// Cast this Quat to a new coordinate type `U`.
    template <typename U>
    inline operator Quat<U>() const {
        return Quat<U>(
            (U)(this->x),
            (U)(this->y),
            (U)(this->z),
            (U)(this->w)
        );
    }
    
    /*******************************
     * Methods                     *
     *******************************/
    
    /// The vector part of this quaternion
    inline Vec<T,3> imag() const {
        return Vec<T,3>(
                detail::VecBase<T,4>::x,
                detail::VecBase<T,4>::y,
                detail::VecBase<T,4>::z);
    }
    
    /// The scalar part of this quaternion
    inline T real() const {
        return detail::VecBase<T,4>::w;
    }
    
    /// The complex conjugate of this quaternion; also the inverse rotation.
    inline Quat<T> conj() const {
        return Quat<T>(-detail::VecBase<T,4>::x, 
                       -detail::VecBase<T,4>::y, 
                       -detail::VecBase<T,4>::z, 
                        detail::VecBase<T,4>::w);
    }
    
    /// Quaternion multiplication
    inline Quat<T> mult(const Quat<T>& q) const {
        Quat<T> result;
        
        const T x = detail::VecBase<T,4>::x;
        const T y = detail::VecBase<T,4>::y;
        const T z = detail::VecBase<T,4>::z;
        const T w = detail::VecBase<T,4>::w;
        
        result.x = w*q.x + x*q.w + y*q.z - z*q.y;
        result.y = w*q.y - x*q.z + y*q.w + z*q.x;
        result.z = w*q.z + x*q.y - y*q.x + z*q.w;
        result.w = w*q.w - x*q.x - y*q.y - z*q.z;
        
        return result;
    }
    
    /** 
     * @brief Convert this unit quaternion to an axis-angle rotation representation
     * `(x, y, z, radians)`.
     */
    Vec<T,4> to_axis_angle_rotation() const {
        T w_clamp = std::min(std::max(detail::VecBase<T,4>::w, -1.0), 1.0);
        double alpha = 2 * std::acos(w_clamp);
        if (detail::VecBase<T,4>::x == 0 and 
            detail::VecBase<T,4>::y == 0 and
            detail::VecBase<T,4>::z == 0)
        {
            return Vec<T,4>(1, 0, 0, 0);
        } else {
            return Vec<T,4>(imag().unit(), alpha);
        }
    }
    
    /**
     * @brief Convert this unit quaternion to an angular velocity representation,
     * with the magnitude of the axis representing the rate of rotation in radians.
     */
    Vec<T,3> to_angular_velocity() const {
        // assumes a unit quaternion.
        T w_clamp = std::min(std::max(detail::VecBase<T,4>::w, -1.0), 1.0);
        T alpha = 2 * std::acos(w_clamp);
        if (imag().is_zero()) {
            return Vec<T,3>::zeros;
        } else {
            return alpha * imag().unit();
        }
    }
    
    /** 
     * @brief Interpolate this unit quaternion with the null rotation according
     * to the interpolation parameter `0 <= t <= 1`.
     */
    // this may go the "long way" around. that what we want?
    Quat<T> rotation_scale(T t) const {
        // assumes a unit quaternion.
        const T w = std::max(-1, std::min(Vec<T,4>::w, 1));
        T omega = std::acos(w);
        T theta = omega * t;
        return Quat<T>(std::sin(theta) * imag().unit(), std::cos(theta));
    }

    /**
     * @brief Spherical linear interpolation, per Ken Shoemake. Interpolate smoothly
     * between two quaternion orientation states.
     * @param q1 The rotation at `t = 1`.
     * @param t Interpolation parameter between 0 and 1.
     * @return 
     */
    Quat<T> slerp(const Quat<T> &q1, T t) const {
        // assumed unit quaternion.
        const Quat<T> &q0 = *this;
        T dot = q0.dot(q1);
        if (dot < 0) { //angle is greater than 90. Avoid "long path" rotations.
            dot = -dot;
            q1 = -q1;
        }
        dot = std::max(0, std::min<T>(1, dot)); //keep <dot> in the domain of acos.
        T omega = std::acos(dot);
        T theta = omega * t;
        Quat<T> q2 = (q1 - q0 * dot).unit();
        return std::cos(theta) * q0 + std::sin(theta) * q2;
    }
    
    /**
     * @brief Quaternion exponential `e`<sup>`q`</sup>.
     * @return A quaternion representing a rotation about `q.imag()` by
     * angle `2 * |q|` and a scaling by `e`<sup>`q.real()`</sup>.
     */
    inline Quat<T> exp() const {
        return std::exp(*this);
    }
    
    /**
     * @brief Quaternion natural log.
     */
    inline Quat<T> log() const {
        return std::log(*this);
    }
    
};


// hashing
template <typename T, typename H>
struct Digest<Quat<T>, H> {
    H operator()(const Quat<T> &v) const {
        H nonce = geom::truncated_constant<H>(0x111113181fd05e73, 0xf5ba45c0011f1ffc);
        return geom::hash_many<H>(nonce, v.x, v.y, v.z, v.w);
    }
};


} //namespace geom

namespace std {
    
/**
 * @brief Quaternion exponential. Represents a rotation about `q.imag()` by
 * angle `|q|` and a scaling by `e`<sup>`q.real()`</sup>.
 */
template <typename T>
geom::Quat<T> exp(const geom::Quat<T>& q) {
    T k = std::exp(q.w);
    T m = q.imag().mag();
    if (m == 0) return geom::Quat<T>(geom::Vec<T,3>::zeros, std::exp(q.w));
    T c = std::cos(m);
    T s = std::sin(m);
    return geom::Quat<T>(k * s * q.imag() / m, k * c);
}

/**
 * @brief Quaternion natural log.
 */
template <typename T>
geom::Quat<T> log(const geom::Quat<T>& q) {
    T q_m = q.mag();
    T v_m = q.imag().mag();
    T k   = std::acos(q.w / q_m);
    return geom::Quat<T>(k * q.imag() / v_m, std::log(q_m));
}

/**
 * @brief Quaternion power. Raise `q` to the (real) power of `a`. 
 */
template <typename T, typename U>
geom::Quat<T> pow(const geom::Quat<T> &q, U a) {
    T mag = q.mag();
    if (mag == 0) return geom::Quat<T>(0,0,0,1);
    U theta = std::acos(q.real() / mag); // :(
    geom::Vec<T,3> nhat = q.imag().unit();
    return std::pow(mag, a) * geom::Quat<T>(
        nhat * std::sin(a * theta), 
        std::cos(a * theta));
}

/**
 * @brief Quaternion complex conjugate. Negate the vector part of `q`. 
 */
template <typename T>
inline geom::Quat<T> conj(const geom::Quat<T> &q) {
    return q.conj();
}

template <typename T>
struct hash<geom::Quat<T>> {
    size_t operator()(const geom::Quat<T> &v) const {
        return geom::hash<geom::Quat<T>, size_t>(v);
    }
};

} // end namespace std

#endif /* QUATERNION_H_ */
