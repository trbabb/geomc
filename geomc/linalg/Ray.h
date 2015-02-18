/*
 * Ray.h
 *
 *  Created on: May 10, 2009
 *      Author: Tim Babb
 */

#ifndef GEOM_RAY_H_
#define GEOM_RAY_H_

#include <geomc/linalg/LinalgTypes.h>
#include <geomc/linalg/Vec.h>

/* RayBase exists so that Ray<T,3> may offer its own functionality
 * but share the rest with the other dimensions.
 */

namespace geom {

namespace detail {

    //ray base class - not for direct use. allows specialization with inheritance
    template <typename T, index_t N> class RayBase {
    public:
        /// Origin of ray.
        Vec<T,N> origin;
        /// Direction of ray.
        Vec<T,N> direction;

        /*************************
         * Structors             *
         *************************/
        
        RayBase():origin((T)0),direction((T)0) {
            direction[0] = 1; //x axis
        }
        
        RayBase(Vec<T,N> o, Vec<T,N> v):origin(o),direction(v) {
            //do nothing else
        }

        /*************************
         * Operators             *
         *************************/

        /// @return A ray pointing in the opposite direction.
        Ray<T,N> operator-() const {
            return neg();
        }

        /*************************
         * Public Methods        *
         *************************/

        /** @return A point along this ray's direction, at distance `d` from
         * the ray origin.
         */
        inline Vec<T,N> atDistance(T d) const {
            return origin + direction.unit()*d;
        }

        /**
         * @return The point `origin + s * direction`.
         */
        inline Vec<T,N> atMultiple(T s) const {
            return origin + direction*s;
        }

        /**
         * @return A ray pointing in the opposite direction.
         */
        inline Ray<T,N> neg() const {
            return Ray<T,N>(origin, -direction);
        }
        
        /** @return The shortest distance from point `p` to the line
         * described by this ray.
         */
        inline T distFromAxisTo(const Vec<T,3> &p) const {
            return std::sqrt(dist2FromAxisTo(p));
        }
        
        /**
         * @return The square of the shortest distance from point `p` to
         * the line described by this ray.
         */
        T dist2FromAxisTo(const Vec<T,3> &p) const {
            Vec<T,N> b = p - origin;
            T b2 = b.mag2();
            T d  = direction.dot(b);
            return b2 - (d * d / b2);
        }

        /**
         * Return the direction normal to `v` which intersects both this ray
         * and the point `p`.
         */
        Vec<T,N> directionFromAxisTo(const Vec<T,N> &p) const {
            Vec<T,N> b = p - origin;
            T d = b.dot(direction);
            return b - (d / direction.mag2()) * direction;
        }
        
    }; //end of RayBase class
    
} // namespace detail

    /** @ingroup linalg
     *  @brief Ray class
     */
    template <typename T, index_t N> 
    class Ray : virtual public detail::RayBase<T,N> {
    public:
        
        /// Construct at the origin, directed down the x+ axis.
        Ray():detail::RayBase<T,N>() {}

        /// Construct a ray with origin `o` and direction `v`.
        Ray(Vec<T,N> o, Vec<T,N> v):detail::RayBase<T,N>(o,v) {}
    };

    /** @ingroup linalg
     *  @brief 3D specialization of Ray class
     */
    template <typename T> class Ray<T,3> : virtual public detail::RayBase<T,3> {
    public:
        /// Construct at the origin, directed down the x+ axis.
        Ray():detail::RayBase<T,3>() {}

        /// Construct a ray with origin `o` and direction `v`.
        Ray(Vec<T,3> o, Vec<T,3> v):detail::RayBase<T,3>(o,v) {}

        // TODO: If we could get a cross product-free implementation of,
        // the below, the entire Ray<3> specialization could be done away with.
        
        /**
         * Return a ray connecting `this` with `r` at their closest approach.
         */
        Vec<T,3> closestApproach(Ray<T,3> r) const {
            Vec<T,3> perp = (r.direction^this->direction).unit();
            return perp * (r.origin - this->origin).dot(perp);
        }

    }; //end of Ray3 specialization

    /**
     * Ray multiple.
     * @return `(r.origin + s * r.direction)`.
     * @related Ray
     */
    template <typename T, index_t N> inline Vec<T,N> operator*(const Ray<T,N> r, T s) {
        return r.atMultiple(s);
    }


    /**
     * Ray multiple.
     * @return `(r.origin + s * r.direction)`.
     * @related Ray
     */
    template <typename T, index_t N> inline Vec<T,N> operator*(T s, const Ray<T,N> r) {
        return r.atMultiple(s);
    }
    
#ifdef GEOMC_LINALG_USE_STREAMS
    
    /** @ingroup linalg
     * Ray stream output, in the form:
     * 
     *     <(x0, x1, x2, ... ), (v0, v1, v2, ...)>
     * where `x` is the origin and `v` is the direction.
     * @related Ray
     */
    template <typename T, index_t N>
    inline std::ostream &operator<< (std::ostream &stream, const Ray<T,N> &r) {
        stream << "<" << r.origin << ", " << r.direction << ">";
        return stream;
    }
    
#endif
} // namespace geom
#endif /* GEOM_RAY_H_ */
