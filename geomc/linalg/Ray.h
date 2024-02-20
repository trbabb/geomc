#ifndef GEOM_RAY_H_
#define GEOM_RAY_H_

#ifdef GEOMC_USE_STREAMS
#include <iostream>
#endif

#include <geomc/linalg/LinalgTypes.h>
#include <geomc/linalg/Vec.h>


namespace geom {

/**
 * @brief Ray class.
 * @tparam T Coordinate type.
 * @tparam N Dimensionality.
 * @ingroup linalg
 */
template <typename T, index_t N>
class Ray {
public:
    /**
     * @brief The type of a point with the same dimension as this ray; a `Vec<T,N>`
     * if `N` is greater than 1; or a `T` if `N` is 1.
     */
    typedef typename PointType<T,N>::point_t point_t;
    /// Origin of ray.
    point_t origin;
    /// Direction of ray.
    point_t direction;

    /*************************
     * Structors             *
     *************************/
    
    /// Construct a ray through the origin directed along the +x axis.
    Ray():origin((T)0),direction((T)0) {
        direction[0] = 1; //x axis
    }
    
    /// Construct a ray through the point `o` with direction `v`.
    Ray(point_t o, point_t v):
        origin(o),
        direction(v) {}

    /*************************
     * Operators             *
     *************************/

    /// Return a ray pointing in the opposite direction.
    Ray<T,N> operator-() const {
        return {origin, -direction};
    }

    /*************************
     * Public Methods        *
     *************************/
    
    /// Change the dimensionality of this Ray.
    template <index_t M>
    inline Ray<T,M> resized() const {
        return Ray<T,M>(
            origin.   template resized<M>(),
            direction.template resized<M>());
    }

    /// Return the point along this ray's direction, at distance `d` from the ray origin.
    inline point_t at_distance(T d) const {
        return origin + direction.unit() * d;
    }

    /// Return the point `origin + s * direction`.
    inline point_t at_multiple(T s) const {
        return origin + direction * s;
    }
    
    /**
     * @brief Find the point on the ray nearest to `p` by orhtogonally
     * projecting `p` onto it.
     */
    inline point_t project(point_t p) const {
        return (p - origin).project_on(direction) + origin;
    }
    
    /// Signed distance from `p` to the nearest point on the ray.
    inline T sdf(point_t p) const {
        return std::sqrt(dist2());
    }
    
    /// Compute the square of the distance from `p` to the nearest point on the ray.
    inline T dist2(point_t p) const {
        point_t b = p - origin;
        T b2 = b.mag2();
        T d  = direction.dot(b);
        return b2 == 0 ? 0 : b2 - (d * d / b2);
    }
    
}; // class Ray


/**
 * Ray multiple.
 * @return `(r.origin + s * r.direction)`.
 * @related Ray
 * @ingroup linalg
 */
template <typename T, index_t N>
inline typename PointType<T,N>::point_t operator*(const Ray<T,N> r, T s) {
    return r.at_multiple(s);
}


/**
 * Ray multiple.
 * @return `(r.origin + s * r.direction)`.
 * @related Ray
 * @ingroup linalg
 */
template <typename T, index_t N>
inline typename PointType<T,N>::point_t operator*(T s, const Ray<T,N> r) {
    return r.at_multiple(s);
}
    
#ifdef GEOMC_USE_STREAMS
    
/** @ingroup linalg
 * Ray stream output, in the form:
 * 
 *     <(x0, x1, x2, ... ), (v0, v1, v2, ...)>
 * where `x` is the origin and `v` is the direction.
 * @related Ray
 */
template <typename T, index_t N>
inline std::ostream &operator<< (std::ostream &stream, const Ray<T,N> &r) {
    stream << "<" << r.origin << " + s * " << r.direction << ">";
    return stream;
}

#endif
} // namespace geom

template <typename T, index_t N>
struct std::hash<geom::Ray<T,N>> {
    size_t operator()(const geom::Ray<T,N> &v) const {
        constexpr size_t nonce = (size_t) 0x8af9e7e642670825ULL;
        return geom::hash_combine(
            geom::hash(v.origin),
            geom::hash(v.direction)
        ) ^ nonce;
    }
};

#endif /* GEOM_RAY_H_ */
