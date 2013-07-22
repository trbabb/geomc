/*
 * Vec2.h
 *
 *  Created on: May 9, 2009
 *      Author: Tim Babb
 */

#ifndef Vec2_H_
#define Vec2_H_

#include <cmath>
#include "VecBase.h"

#ifdef GEOMC_VEC_USE_SWIZZLE
    #include <string>
    #include "GeomException.h"
#endif

namespace geom {

    //Partial specialization of Vec<T,N> for N = 2
    template <typename T> class Vec<T,2> : public detail::VecCommon<T,2> {

    public:
        
        static const Vec<T,2> X_AXIS;
        static const Vec<T,2> Y_AXIS;

        /*===========================*
         * Structors                 *
         *===========================*/
        
        Vec():detail::VecCommon<T,2>(){}
        
        Vec(T a):detail::VecCommon<T,2>(a){}

        Vec(T x, T y){
            detail::VecBase<T,2>::x = x;
            detail::VecBase<T,2>::y = y;
        }
        
        Vec(T v[2]):detail::VecCommon<T,2>(v){}
        
        /*===========================*
         * Convenience Arithmetic    *
         *===========================*/
        
        //do not hide these; we want to overload them:
        using detail::VecCommon<T,2>::add;
        using detail::VecCommon<T,2>::minus;
        using detail::VecCommon<T,2>::scale;
        using detail::VecCommon<T,2>::dot;

        inline Vec<T,2> add(T dx, T dy) const {
            return Vec<T,2>(detail::VecBase<T,2>::x + dx, detail::VecBase<T,2>::y + dy);
        }

        inline Vec<T,2> minus(T dx, T dy) const {
            return Vec<T,2>(detail::VecBase<T,2>::x - dx, detail::VecBase<T,2>::y - dy);
        }

        inline Vec<T,2> scale(T sx, T sy) const {
            return Vec<T,2>(detail::VecBase<T,2>::x*sx, detail::VecBase<T,2>::y*sy);
        }

        inline T dot(T x1, T y1) const {
            return detail::VecBase<T,2>::x*x1 + detail::VecBase<T,2>::y*y1;
        }

        inline Vec<T,2> reflectX() const {
            return Vec<T,2>(-detail::VecBase<T,2>::x, detail::VecBase<T,2>::y);
        }

        inline Vec<T,2> reflectY() const {
            return Vec<T,2>(detail::VecBase<T,2>::x, -detail::VecBase<T,2>::y);
        }

        bool isMultipleOf(Vec<T,2> v) const {
            T ratio = detail::VecBase<T,2>::x/v.x;
            return ratio == (detail::VecBase<T,2>::y/v.y);
        }

        bool isFiniteReal() const {
            return isReal(detail::VecBase<T,2>::x) && isReal(detail::VecBase<T,2>::y);
        }
        
        /*===========================*
         * 2D Geometry               *
         *===========================*/

        Vec<T,2> leftPerpendicular() const {
            return Vec<T,2>(-detail::VecBase<T,2>::y, detail::VecBase<T,2>::x);
        }

        Vec<T,2> rightPerpendicular() const {
            return Vec<T,2>(detail::VecBase<T,2>::y, -detail::VecBase<T,2>::x);
        }

        Vec<T,2> rotate(T radians) const {
            T sint = sin(radians);
            T cost = cos(radians);
            return Vec<T,2>(cost*detail::VecBase<T,2>::x - sint*detail::VecBase<T,2>::y,
                            sint*detail::VecBase<T,2>::x + cost*detail::VecBase<T,2>::y);
        }

        Vec<T,2> rotate(Vec<T,2> center, T radians) const {
            return center + (((*this) - center).rotate(radians));
        }

        //returns (r, theta) from (x, y)
        inline Vec<T,2> toPolar() const {
            return Vec<T,2>(this->mag(), getAngle());
        }
        
        //returns (x, y) from (r, theta)
        inline Vec<T,2> fromPolar() const {
            Vec<T,2> v;
            v.x = detail::VecBase<T,2>::y * cos(detail::VecBase<T,2>::x);
            v.y = detail::VecBase<T,2>::y * sin(detail::VecBase<T,2>::x);
            return v;
        }

        inline T getAngle() const {
            T theta = (T)(atan2(detail::VecBase<T,2>::y,detail::VecBase<T,2>::x));
            if (theta < 0) theta += M_2_PI;
            return theta;
        }
        
        #ifdef GEOMC_VEC_USE_SWIZZLE
        const Vec<T,2> swizzle(const std::string& xyz) const {
            if (xyz.length() != 2){
                throw GeomException("Invalid 2D swizzle string: " + xyz);
            }
            return Vec<T,2>(swizzlecoord(xyz[0]),
                            swizzlecoord(xyz[1]));
        }
        #endif

    private:
        
        #ifdef GEOMC_VEC_USE_SWIZZLE
        T swizzlecoord(char c) const  {
            switch (tolower(c)){
                case 'x':
                    return detail::VecBase<T,2>::x;
                case 'y':
                    return detail::VecBase<T,2>::y;
                case '0':
                    return 0;
                case '1':
                    return 1;
                default:
                    throw GeomException("Invalid swizzle character: " + c);
            }
        }
        #endif

    }; //end Vec2 definition
    
    template <typename T> const Vec<T,2> Vec<T,2>::X_AXIS = Vec<T,2>(1,0);
    template <typename T> const Vec<T,2> Vec<T,2>::Y_AXIS = Vec<T,2>(0,1);

}; //end namespace geom

#endif /* Vec2_H_ */
