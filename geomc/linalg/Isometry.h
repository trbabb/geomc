#pragma once

#include <geomc/linalg/Rotation.h>

namespace geom {

/**
 * @ingroup linalg
 * @brief A rigid rotation and translation.
 * 
 * Isometric transfroms do not have any skew or scales; they preserve
 * shapes, angles, and distances.
 */
template <typename T, index_t N>
struct Isometry {
    /// Rotation component.
    Rotation<T,N> rx;
    /// Translation component.
    Vec<T,N>      tx;
    
    Isometry():rx(),tx() {}
    Isometry(const Rotation<T,N>& rx, const Vec<T,N>& tx):rx(rx),tx(tx) {}
    Isometry(const Rotation<T,N>& rx):rx(rx),tx() {}
    Isometry(const Vec<T,N>& tx):rx(),tx(tx) {}
    
    /// Cast to an affine transform.
    operator AffineTransform<T,N>() const {
        return geom::translation(tx) * rx.transform();
    }
    
    /// Compose two isometric transforms.
    Isometry<T,N> operator*(const Isometry<T,N>& other) const {
        return Isometry<T,N>(rx * other.rx, tx + rx.transform(other.tx));
    }
    
    /// Transform a point.
    Vec<T,N> operator*(const Vec<T,N>& p) const {
        return rx * p + tx;
    }
    
    /// Transform a direction vector.
    Vec<T,N> apply_direction(const Vec<T,N> v) const {
        return rx * v;
    }
    
    /// Inverse transform a direction vector.
    Vec<T,N> apply_inverse_direction(const Vec<T,N> v) const {
        return v / rx;
    }
    
    /// Compose with the inverse of an isometry.
    Isometry<T,N> operator/(const Isometry<T,N>& other) const {
        return Isometry<T,N>(rx / other.rx, rx.transform(other.tx - tx));
    }
    
    /// Compute the inverse of the isometry.
    Isometry<T,N> inverse() const {
        Rotation<T,N> r_inv = rx.inverse();
        return Isometry<T,N>(
            r_inv,
            r_inv * (-tx)
        );
    }
    
    /// Apply a translation.
    Isometry<T,N> operator+=(const Vec<T,N>& v) {
        tx += v;
        return *this;
    }
    
    /// Compose in-place.
    Isometry<T,N> operator*=(const Isometry<T,N>& other) {
        *this = *this * other;
        return *this;
    }
};

/** @addtogroup linalg
 *  @{
 */

/// @brief Transform a point.
/// @related Isometry
template <typename T, index_t N>
Vec<T,N> operator/(const Vec<T,N>& v, const Isometry<T,N>& i) {
    // unapply the translation and rotation in reverse order
    return (v - i.tx) / i.rx;
}

/// @brief Transform a ray.
/// @related Isometry
/// @related Ray
template <typename T, index_t N>
Ray<T,N> operator*(const Isometry<T,N>& xf, const Ray<T,N>& ray) {
    return {xf * ray.origin, xf.rx * ray.direction};
}

/// @brief Inverse-transform a ray.
/// @related Isometry
/// @related Ray
template <typename T, index_t N>
Ray<T,N> operator/(const Ray<T,N>& ray, const Isometry<T,N>& xf) {
    return {ray.origin / xf, ray.direction / xf.rx};
}

/// @brief Apply a rotation to an isometry.
/// @related Isometry
/// @related Rotation
template <typename T, index_t N>
Isometry<T,N> operator*(const Rotation<T,N>& r, const Isometry<T,N>& i) {
    return Isometry<T,N>(r * i.rx, r.transform(i.tx));
}

/// @brief Apply an isometry to a rotation.
/// @related Isometry
/// @related Rotation
template <typename T, index_t N>
Isometry<T,N> operator*(const Isometry<T,N>& i, const Rotation<T,N>& r) {
    return Isometry<T,N>(i.rx * r, i.tx);
}

/// @brief Apply a translation to an isometry.
/// @related Isometry
template <typename T, index_t N>
Isometry<T,N> operator+(const Isometry<T,N>& i, const Vec<T,N>& v) {
    return Isometry<T,N>(i.rx, i.tx + v);
}

/// @brief Apply a translation to an isometry.
/// @related Isometry
template <typename T, index_t N>
Isometry<T,N> operator+(const Vec<T,N>& v, const Isometry<T,N>& i) {
    return Isometry<T,N>(i.rx, i.tx + v);
}

/// @brief Scale the magnitude of an isometry. A scale of 0 produces an identity transform.
/// Applying a scale of 1 to an isometry results in no change.
/// @related Isometry
template <typename T, index_t N>
Isometry<T,N> operator*(T s, const Isometry<T,N>& xf) {
    return {
        s * xf.rx,
        s * xf.tx
    };
}

/// @brief Scale the magnitude of an isometry.
/// @related Isometry
template <typename T, index_t N>
Isometry<T,N> operator*(const Isometry<T,N>& xf, T s) {
    return {
        s * xf.rx,
        s * xf.tx
    };
}

/// @brief Continuously interpolate two isometries.
/// @related Isometry
template <typename T, index_t N>
Isometry<T,N> mix(T s, const Isometry<T,N>& a, const Isometry<T,N>& b) {
    return {
        mix(s, a.rx, b.rx),
        mix(s, a.tx, b.tx)
    };
}

/// @}  // group linalg

template <typename T, index_t N, typename H>
struct Digest<Isometry<T,N>, H> {
    H operator()(const Isometry<T,N>& xf) const {
        H nonce = geom::truncated_constant<H>(0x176a7504edd40424, 0x292b420bad0c4b61);
        return geom::hash_many(nonce, xf.rx, xf.tx);
    }
};

} // namespace geom
