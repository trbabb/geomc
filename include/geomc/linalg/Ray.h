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
        Vec<T,N> origin;
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

        Ray<T,N> operator-() const {
            return neg();
        }

        /*************************
         * Public Methods        *
         *************************/

        inline Vec<T,N> atDistance(T d) const {
            return origin + direction.unit()*d;
        }

        inline Vec<T,N> atMultiple(T s) const {
            return origin + direction*s;
        }

        inline Ray<T,N> neg() const {
            return Ray<T,N>(origin, -direction);
        }
    }; //end of RayBase class
    
}; // namespace detail

    //Ray class
    template <typename T, index_t N> class Ray : virtual public detail::RayBase<T,N> {
    public:
        Ray():detail::RayBase<T,N>() {}

        Ray(Vec<T,N> o, Vec<T,N> v):detail::RayBase<T,N>(o,v) {}
    };

    //Ray3 specialization
    template <typename T> class Ray<T,3> : virtual public detail::RayBase<T,3> {
    public:
        Ray():detail::RayBase<T,3>() {}

        Ray(Vec<T,3> o, Vec<T,3> v):detail::RayBase<T,3>(o,v) {}

        T distFromAxisTo(Vec<T,3> v) const {
            return std::sqrt((this->direction^(this->origin-v)).mag2() / this->direction.mag2());
        }

        Vec<T,3> directionFromAxisTo(Vec<T,3> v) const {
            Vec<T,3> perp = this->direction^(v^this->direction);
            return perp * (v.dot(perp)/perp.mag());
        }

        Vec<T,3> closestApproach(Ray<T,3> r) const {
            Vec<T,3> perp = (r.direction^this->direction).unit();
            return perp * (r.origin - this->origin).dot(perp);
        }

    }; //end of Ray3 specialization


    template <typename T, index_t N> inline Vec<T,N> operator*(const Ray<T,N> r, T s) {
        return r.atMultiple(s);
    }


    template <typename T, index_t N> inline Vec<T,N> operator*(T s, const Ray<T,N> r) {
        return r.atMultiple(s);
    }
    
#ifdef GEOMC_LINALG_USE_STREAMS
    
    template <typename T, index_t N>
    inline std::ostream &operator<< (std::ostream &stream, const Ray<T,N> &r) {
        stream << "<" << r.origin << ", " << r.direction << ">";
        return stream;
    }
    
#endif
};
#endif /* GEOM_RAY_H_ */
