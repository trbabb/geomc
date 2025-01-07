#pragma once

#include <geomc/linalg/Rotation.h>

// todo: lie-algebra operations (exp, log, etc)

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
    Isometry<T,N> operator+=(const Vec<T,N>& v) {
        tx += v;
        return *this;
    }
    
    Isometry<T,N> operator*=(const Isometry<T,N>& other) {
        *this = *this * other;
        return *this;
    }
};

template <typename T, index_t N>
Isometry<T,N> operator*(const Rotation<T,N>& r, const Isometry<T,N>& i) {
    return Isometry<T,N>(r * i.rx, r.transform(i.tx));
}

template <typename T, index_t N>
Isometry<T,N> operator*(const Isometry<T,N>& i, const Rotation<T,N>& r) {
    return Isometry<T,N>(i.rx * r, i.tx);
}

template <typename T, index_t N>
Isometry<T,N> operator+(const Isometry<T,N>& i, const Vec<T,N>& v) {
    return Isometry<T,N>(i.rx, i.tx + v);
}

template <typename T, index_t N>
Isometry<T,N> operator+(const Vec<T,N>& v, const Isometry<T,N>& i) {
    return Isometry<T,N>(i.rx, i.tx + v);
}


} // namespace geom
