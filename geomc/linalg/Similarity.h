#pragma once

#include <geomc/linalg/Isometry.h>

// todo: lie algebra operations (exp, log, etc)

namespace geom {

/**
 * @ingroup linalg
 * @brief A similarity transform, which is a rotation, scaling, and translation.
 *
 * Similar transfroms do not have any skew or nonuniform scales; they preserve
 * shapes, angles, and relative distances.
 *
 * For transforms which do not include a scaling, see Isometry.
 *
 * For nonuniform scaling or skew transforms, see AffineTransform.
 *
 * Similarity transforms compose like transforms, with multiplication on the left:
 *
 *     Similarity<T,N> s1, s2;
 *     Similarity<T,N> s3 = s2 * s1; // similarity which applies s2, then s1
 *     Vec<T,N> v = s3 * s1 * v0;    // apply s1 to v0, then s3 to the result
 *     Sphere<T,N> s = s1 * sphere;  // apply s1 to a sphere
 *
 * Similarity transforms can be inverted with the `/` operator:
 *
 *     Similarity<T,N> s1, s2;
 *     Similarity<T,N> s3 = s2 / s1; // similarity which takes s1 to s2
 *     Vec<T,N> v = v0 / s3;       // apply the inverse of s3 to v0
 *
 * Compose with rotations and translations:
 *
 *     Similarity<T,3> s;
 *     Rotation<T,3> r;
 *     Similarity<T,3> s2 = r * s; // s2 is s with r applied
 *     Similarity<T,3> s3 = s + Vec<T,3>(1,2,3); // s3 is s translated by (1,2,3)
 *
 * Scaling is done explicitly with `.scaled()` or the `.sx` member:
 *
 *     Similarity<T,3> s;
 *     Similarity<T,3> s2 = s.scaled(2); // s2 is s with double the scale
 *     Similarity<T,3> s3 = s * 2; // s3 is s with double the scale
 */
template <typename T, index_t N>
struct Similarity {
    /// Scale component.
    T             sx = 1;
    /// Rotation component.
    Rotation<T,N> rx;
    /// Translation component.
    Vec<T,N>      tx;
    
    Similarity():rx(),tx() {}
    Similarity(T sx, const Rotation<T,N>& rx, const Vec<T,N>& tx):sx(sx),rx(rx),tx(tx) {}
    Similarity(T sx, const Rotation<T,N>& rx):sx(sx),rx(rx),tx() {}
    Similarity(const Rotation<T,N>& rx):rx(rx),tx() {}
    Similarity(const Vec<T,N>& tx):rx(),tx(tx) {}
    Similarity(const Isometry<T,N>& s, T sx=1):sx(sx),rx(s.rx),tx(s.tx) {}
    explicit Similarity(T sx):sx(sx),rx(),tx() {}
    
    /// Represent this similarity as an affine transform.
    operator AffineTransform<T,N>() const {
        return geom::translation(tx) * rx.transform() * geom::scale(sx);
    }
    
    /// Cast the underlying coordinate type.
    template <typename U>
    explicit operator Similarity<U,N>() const {
        return Similarity<U,N>(sx, rx, tx);
    }
    
    /// Extend the dimensionality of this similarity.
    template <index_t M>
    requires (M > N)
    operator Similarity<T,M>() const {
        return Similarity<T,M>(sx, rx, tx.template resized<M>());
    }
    
    /// Extend the dimensionality of this similarity and change the coordinate type.
    template <typename U, index_t M>
    requires (M > N)
    explicit operator Similarity<U,M>() const {
        return Similarity<U,M>(
            static_cast<Rotation<U,M>>(rx),
            Vec<U,M>(tx)
        );
    }
    
    /// Transform a point.
    Vec<T,N> operator*(const Vec<T,N>& p) const {
        return tx + sx * rx * p;
    }
    
    /// Transform a direction vector.
    Vec<T,N> apply_direction(const Vec<T,N>& v) const {
        return sx * rx.transform(v);
    }
    
    /// Inverse-transform a direction vector.
    Vec<T,N> apply_inverse_direction(const Vec<T,N>& v) const {
        return (v / sx) / rx;
    }
    
    /// Compose two similarity transforms.
    Similarity<T,N> operator*(const Similarity<T,N>& other) const {
        return Similarity<T,N>(
            sx * other.sx,
            rx * other.rx,
            tx + rx.transform(other.tx) * sx
        );
    }
    
    /// Compose with the inverse of a similarity.
    Similarity<T,N> operator/(const Similarity<T,N>& other) const {
        return Similarity<T,N>(
            sx / other.sx,
            rx / other.rx,
            rx.transform(other.tx - tx) / other.sx
        );
    }
    
    /// Compute the inverse of the similarity.
    Similarity<T,N> inverse() const {
        Rotation<T,N> r_inv = rx.inverse();
        return Similarity<T,N>(
            1 / sx,
            r_inv,
            r_inv * (-tx) / sx
        );
    }
    
    Similarity<T,N> operator+=(const Vec<T,N>& v) {
        tx += v;
        return *this;
    }
    
    Similarity<T,N> operator*=(const Similarity<T,N>& other) {
        *this = *this * other;
        return *this;
    }
    
    Similarity<T,N> scaled(T s) const {
        return Similarity<T,N>(s * sx, rx, tx);
    }
    
};

/** @addtogroup linalg
 * @{
 */

/// @brief Transform a point.
/// @related Similarity
template <typename T, index_t N>
Vec<T,N> operator/(const Vec<T,N>& v, const Similarity<T,N>& i) {
    // unapply the translation, rotation, and scale in reverse order
    return ((v - i.tx) / i.rx) / i.sx;
}

/// @brief Transform a ray.
/// @related Similarity
/// @related Ray
template <typename T, index_t N>
Ray<T,N> operator*(const Similarity<T,N>& i, const Ray<T,N>& ray) {
    return {i * ray.origin, i.rx * i.sx * ray.direction};
}

/// @brief Inverse-transform a ray.
/// @related Similarity
/// @related Ray
template <typename T, index_t N>
Ray<T,N> operator/(const Ray<T,N>& ray, const Similarity<T,N>& i) {
    return {ray.origin / i, (ray.direction / i.sx) / i.rx};
}

/// @brief Apply a rotation to a similarity.
/// @related Similarity
/// @related Rotation
template <typename T, index_t N>
Similarity<T,N> operator*(const Rotation<T,N>& r, const Similarity<T,N>& i) {
    return Similarity<T,N>(i.sx, r * i.rx, r.transform(i.tx));
}

/// @brief Apply a similarity to a rotation.
/// @related Similarity
/// @related Rotation
template <typename T, index_t N>
Similarity<T,N> operator*(const Similarity<T,N>& i, const Rotation<T,N>& r) {
    return Similarity<T,N>(i.sx, i.rx * r, i.tx);
}

/// @brief Apply a translation to a similarity.
/// @related Similarity
template <typename T, index_t N>
Similarity<T,N> operator+(const Vec<T,N>& v, const Similarity<T,N>& i) {
    return Similarity<T,N>(i.sx, i.rx, i.tx + v);
}

/// @brief Apply a translation to a similarity.
/// @related Similarity
template <typename T, index_t N>
Similarity<T,N> operator+(const Similarity<T,N>& i, const Vec<T,N>& v) {
    return Similarity<T,N>(i.sx, i.rx, i.tx + v);
}

/// @brief Interpolate between two similarity transforms.
/// It is invalid to interpolate between two similarity transforms with different signs.
template <typename T, index_t N>
Similarity<T,N> mix(T s, const Similarity<T,N>& a, const Similarity<T,N>& b) {
    T sign = a.sx >= 0 ? 1 : -1;
    return Similarity<T,N>(
        sign * std::pow(std::abs(a.sx), 1 - s) * std::pow(std::abs(b.sx), s),
        geom::mix(s, a.rx, b.rx),
        geom::mix(s, a.tx, b.tx)
    );
}

// nb: we omit operator* with a scalar for now; it's ambiguous whether
//   it means to apply a scaling or a partial application.


/// @} // end addtogroup linalg

template <typename T, index_t N, typename H>
struct Digest<Similarity<T,N>, H> {
    H operator()(const Isometry<T,N>& xf) const {
        H nonce = geom::truncated_constant<H>(0xcab7aaee307c0de5, 0xd3074f5dbd666d01);
        return geom::hash_many(nonce, xf.rx, xf.tx, xf.sx);
    }
};

} // namespace geom
