#pragma once

#include <geomc/linalg/Rotation.h>

// todo: lie algebra operations (exp, log, etc)

namespace geom {

template <typename T, index_t N>
struct Similarity {
    T             sx = 1;
    Rotation<T,N> rx;
    Vec<T,N>      tx;
    
    Similarity():rx(),tx() {}
    Similarity(T sx, const Rotation<T,N>& rx, const Vec<T,N>& tx):sx(sx),rx(rx),tx(tx) {}
    Similarity(T sx, const Rotation<T,N>& rx):sx(sx),rx(rx),tx() {}
    Similarity(const Rotation<T,N>& rx):rx(rx),tx() {}
    Similarity(const Vec<T,N>& tx):rx(),tx(tx) {}
    explicit Similarity(T sx):sx(sx),rx(),tx() {}
    
    /// Represent this similarity as an affine transform.
    operator AffineTransform<T,N>() const {
        return geom::translation(tx) * rx.transform() * geom::scale(sx);
    }
    
    /// Apply the similarity to a point.
    Vec<T,N> operator*(const Vec<T,N>& p) const {
        return tx + sx * rx.transform(p);
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

template <typename T, index_t N>
Similarity<T,N> operator*(const Rotation<T,N>& r, const Similarity<T,N>& i) {
    return Similarity<T,N>(i.sx, r * i.rx, r.transform(i.tx));
}

template <typename T, index_t N>
Similarity<T,N> operator*(const Similarity<T,N>& i, const Rotation<T,N>& r) {
    return Similarity<T,N>(i.sx, i.rx * r, i.tx);
}

template <typename T, index_t N>
Similarity<T,N> operator+(const Vec<T,N>& v, const Similarity<T,N>& i) {
    return Similarity<T,N>(i.sx, i.rx, i.tx + v);
}

template <typename T, index_t N>
Similarity<T,N> operator+(const Similarity<T,N>& i, const Vec<T,N>& v) {
    return Similarity<T,N>(i.sx, i.rx, i.tx + v);
}

// nb: we omit operator* denoting scaling, because it might be confused with
//    a partial or iterated application of the similarity transform.


} // namespace geom
