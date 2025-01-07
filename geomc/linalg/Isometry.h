#pragma once

#include <geomc/linalg/Rotation.h>

namespace geom {

/**
 * A rigid rotation and translation.
 */
template <typename T, index_t N>
struct Isometry {
    Rotation<T,N> rx;
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
    
    /// Apply the isometry to a point.
    Vec<T,N> operator*(const Vec<T,N>& p) const {
        return rx.transform(p) + tx;
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

/// Apply a rotation to an isometry.
template <typename T, index_t N>
Isometry<T,N> operator*(const Rotation<T,N>& r, const Isometry<T,N>& i) {
    return Isometry<T,N>(r * i.rx, r.transform(i.tx));
}

/// Apply an isometry to a rotation.
template <typename T, index_t N>
Isometry<T,N> operator*(const Isometry<T,N>& i, const Rotation<T,N>& r) {
    return Isometry<T,N>(i.rx * r, i.tx);
}

/// Apply a translation to an isometry.
template <typename T, index_t N>
Isometry<T,N> operator+(const Isometry<T,N>& i, const Vec<T,N>& v) {
    return Isometry<T,N>(i.rx, i.tx + v);
}

/// Apply a translation to an isometry.
template <typename T, index_t N>
Isometry<T,N> operator+(const Vec<T,N>& v, const Isometry<T,N>& i) {
    return Isometry<T,N>(i.rx, i.tx + v);
}

/// Scale the magnitude of an isometry. A scale of 0 produces an identity transform.
/// Applying a scale of 1 to an isometry results in no change.
template <typename T, index_t N>
Isometry<T,N> operator*(T s, const Isometry<T,N>& xf) {
    return {
        s * xf.rx,
        s * xf.tx
    };
}

/// Scale the magnitude of an isometry.
template <typename T, index_t N>
Isometry<T,N> operator*(const Isometry<T,N>& xf, T s) {
    return {
        s * xf.rx,
        s * xf.tx
    };
}

/// Continuously interpolate two isometries.
template <typename T, index_t N>
Isometry<T,N> mix(T s, const Isometry<T,N>& a, const Isometry<T,N>& b) {
    return {
        mix(s, a.rx, b.rx),
        mix(s, a.tx, b.tx)
    };
}

} // namespace geom
