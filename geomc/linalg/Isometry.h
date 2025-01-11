#pragma once

#include <geomc/linalg/Rotation.h>

namespace geom {

/**
 * @ingroup linalg
 * @brief A rigid rotation and translation.
 * 
 * Isometric transfroms do not have any skew or scales; they preserve
 * shapes, angles, and distances.
 *
 * Isometries meet the Transform concept. For transforms which include a scaling,
 * see Similarity. For nonuniform scaling or skew transforms, see AffineTransform.
 *
 * Isometries compose like transforms, with multiplication on the left:
 *
 *     Isometry<T,N> i1, i2;
 *     Isometry<T,N> i3 = i2 * i1;  // isometry which applies i2, then i1
 *     Vec<T,N> v = i3 * i1 * v0;   // apply i1 to v0, then i3 to the result
 *     Sphere<T,N> s = i1 * sphere; // apply i1 to a sphere
 *
 * Isometries can be inverted with the `/` operator:
 *
 *     Isometry<T,N> i1, i2;
 *     Isometry<T,N> i3 = i2 / i1; // isometry which takes i1 to i2
 *     Vec<T,N> v = v0 / i3;       // apply the inverse of i3 to v0
 *
 * Compose with rotations and translations:
 * 
 *     Isometry<T,3> i;
 *     Rotation<T,3> r;
 *     Isoemtry<T,3> i2 = r * i; // i2 is i with r post-applied
 *     Isometry<T,3> i3 = i + Vec<T,3>(1,2,3); // i3 is i translated by (1,2,3)
 * 
 */
template <typename T, index_t N>
class Isometry : public Dimensional<T,N> {
public:
    
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
    
    /// Cast the underlying coordinate type.
    template <typename U>
    explicit operator Isometry<U,N>() const {
        return Isometry<U,N>(rx, tx);
    }
    
    /// Extend the dimensionality of this isometry.
    template <index_t M>
    requires (M > N)
    operator Isometry<T,M>() const {
        return Isometry<T,M>(rx, tx.template resized<M>());
    }
    
    /// Extend the dimensionality of this isometry and change the coordinate type.
    template <typename U, index_t M>
    requires (M > N)
    explicit operator Isometry<U,M>() const {
        return Isometry<U,M>(
            static_cast<Rotation<U,M>>(rx),
            Vec<U,M>(tx)
        );
    }
    
    /// Compose two isometric transforms.
    Isometry<T,N> operator*(const Isometry<T,N>& other) const {
        return Isometry<T,N>(rx * other.rx, tx + rx.transform(other.tx));
    }
    
    /// Compose in-place.
    Isometry<T,N>& operator*=(const Isometry<T,N>& other) {
        tx += rx.transform(other.tx);
        rx *= other.rx;
        return *this;
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
    
    /// In-place apply inverse
    Isometry<T,N>& operator/=(const Isometry<T,N>& other) {
        tx = rx.transform(other.tx - tx);
        rx /= other.rx;
        return *this;
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
