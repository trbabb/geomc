/*
 * Vec3.h
 *
 *  Created on: Feb 22, 2009
 *      Author: Tim Babb
 */

#ifndef VECTOR3D_H_
#define VECTOR3D_H_

#include <cctype>
#include "VecBase.h"

#ifdef GEOMC_VEC_USE_SWIZZLE
    #include <string>
    #include "GeomException.h"
#endif


namespace geom {

    /*===========================*
     * Global Functions          *
     *===========================*/

    template <typename T> 
    inline bool isReal(T d){
        return !(std::isinf(d) || std::isnan(d));
    }

    template <typename T> 
    Vec<T,3> triangle_normal(const Vec<T,3> v1, const Vec<T,3> v2, const Vec<T,3> v3){
        return (v2-v1).cross(v3-v1).unit();
    }

    template <typename T> 
    Vec<T,3> rainbow(double s){
        s = std::fmod(s,1);
        if (s < 0) s = 1+s;
        return Vec<T,3>(
                std::min( std::max( std::max( -6*s+2, 6*s-4 ), 0.0), 1.0),
                std::max( std::min( std::min( 6*s,   -6*s+4 ), 1.0), 0.0),
                std::max( std::min( std::min( 6*s-2, -6*s+6 ), 1.0), 0.0));
    }

    template <typename T> 
    Vec<T,3> fromRGB(int aRGB){
        T r = ((aRGB & 0x00ff0000) >> 16) / detail::_RGBChannelConversion<T>::factor;
        T g = ((aRGB & 0x0000ff00) >> 8)  / detail::_RGBChannelConversion<T>::factor;
        T b =  (aRGB & 0x000000ff)        / detail::_RGBChannelConversion<T>::factor;
        return Vec<T,3>(r,g,b);
    }

    /*===========================*
     * Vector Class              *
     *===========================*/

    template <typename T> class Vec<T,3> : public detail::VecCommon<T,3> {
    public:
        
        static const Vec<T,3> X_AXIS;
        static const Vec<T,3> Y_AXIS;
        static const Vec<T,3> Z_AXIS;

        /*===========================*
         * Structors                 *
         *===========================*/

        Vec():detail::VecCommon<T,3>(){}
        
        //Vec(Vec<T,2> xy, T z):detail::VecBase<T,3>::x(xy.x),detail::VecBase<T,3>::y(xy.y),detail::VecBase<T,3>::z(z){}
        Vec(Vec<T,2> xy, T z){
            detail::VecBase<T,3>::x = xy.x;
            detail::VecBase<T,3>::y = xy.y;
            detail::VecBase<T,3>::z = z;
        }

        Vec(T a):detail::VecCommon<T,3>(a){}

        Vec(T v[3]):detail::VecCommon<T,3>(v){}

        //Vec(T x, T y, T z):detail::VecBase<T,3>::x(x),detail::VecBase<T,3>::y(y),detail::VecBase<T,3>::z(z){}
        Vec(T x, T y, T z){
            detail::VecBase<T,3>::x = x;
            detail::VecBase<T,3>::y = y;
            detail::VecBase<T,3>::z = z;
        }

        /*===========================*
         * Operators                 *
         *===========================*/

        //cross product
        const Vec<T,3>  operator^(Vec<T,3> v) const {
            return cross(v);
        }
        
        /*===========================*
         * Arithmetic                *
         *===========================*/
        
        //do not hide these; we want to overload them: 
        using detail::VecCommon<T,3>::add;
        using detail::VecCommon<T,3>::minus;
        using detail::VecCommon<T,3>::scale;
        using detail::VecCommon<T,3>::dot;
        
        const Vec<T,3> add(T dx, T dy, T dz) const {
            return Vec<T,3>(detail::VecBase<T,3>::x+dx, detail::VecBase<T,3>::y+dy, detail::VecBase<T,3>::z+dz);
        }

        const Vec<T,3> minus(T dx, T dy, T dz) const {
            return Vec<T,3>(detail::VecBase<T,3>::x-dx, detail::VecBase<T,3>::y-dy, detail::VecBase<T,3>::z-dz);
        }

        const Vec<T,3> scale(T a, T b, T c) const {
            return Vec<T,3>(a*detail::VecBase<T,3>::x, b*detail::VecBase<T,3>::y, c*detail::VecBase<T,3>::z);
        }
        
        T dot(T xx, T yy, T zz) const {
            return (detail::VecBase<T,3>::x * xx) + (detail::VecBase<T,3>::y*yy) + (detail::VecBase<T,3>::z*zz);
        }

        const Vec<T,3> reflectX() const {
            return Vec<T,3>(-detail::VecBase<T,3>::x, detail::VecBase<T,3>::y, detail::VecBase<T,3>::z);
        }

        const Vec<T,3> reflectY() const {
            return Vec<T,3>(detail::VecBase<T,3>::x, -detail::VecBase<T,3>::y, detail::VecBase<T,3>::z);
        }

        const Vec<T,3> reflectZ() const {
            return Vec<T,3>(detail::VecBase<T,3>::x, detail::VecBase<T,3>::y, -detail::VecBase<T,3>::z);
        }
        
        /*===========================*
         * 3D Geometry               *
         *===========================*/
        
        const Vec<T,3> cross(Vec<T,3> v) const {
            T newX = (detail::VecBase<T,3>::y * v.z) - (detail::VecBase<T,3>::z * v.y);
            T newY = (detail::VecBase<T,3>::z * v.x) - (detail::VecBase<T,3>::x * v.z);
            T newZ = (detail::VecBase<T,3>::x * v.y) - (detail::VecBase<T,3>::y * v.x);

            return Vec<T,3>(newX, newY, newZ);
        }
        
        //returns radius, elevation, azimuth
        //elevation will be in the range (0, pi)
        //azimuth will be in the range (-pi, pi)
        const Vec<T,3> toSpherical() const {
            T r = this->mag();
            return Vec<T,3>(r, acos(detail::VecBase<T,3>::z/r), atan2(detail::VecBase<T,3>::y,detail::VecBase<T,3>::x));
        }
        
        //returns x,y,z from radius, elevation, azimuth
        const Vec<T,3> fromSpherical() const {
            Vec<T,3> oot;
            T rsinelev = detail::VecBase<T,3>::x * sin(detail::VecBase<T,3>::y);
            oot.x = rsinelev*cos(detail::VecBase<T,3>::z);
            oot.y = rsinelev*sin(detail::VecBase<T,3>::z);
            oot.z = detail::VecBase<T,3>::x*cos(detail::VecBase<T,3>::y);
            return oot;
        }

        //!!!const Vec<T,3> rotate(Vector4d quat) const;

        const Vec<T,3> rotate(Vec<T,3> axis, T radians) const {
            return rotate(axis.x, axis.y, axis.z, radians);
        }

        const Vec<T,3> rotate(T u, T v, T w, T radians) const {
            T u2 = u*u;
            T v2 = v*v;
            T w2 = w*w;
            T uxvywz = u*detail::VecBase<T,3>::x+v*detail::VecBase<T,3>::y+w*detail::VecBase<T,3>::z;
            T uvw2 = u2+v2+w2;
            T rsin = sqrt(uvw2)*sin(radians);
            T c    = cos(radians);

            return Vec<T,3>(
                    (u*uxvywz+(detail::VecBase<T,3>::x*(v2 + w2)-u*(v*detail::VecBase<T,3>::y+w*detail::VecBase<T,3>::z))*c+rsin*(-w*detail::VecBase<T,3>::y+v*detail::VecBase<T,3>::z))/uvw2,
                    (v*uxvywz+(detail::VecBase<T,3>::y*(u2 + w2)-v*(u*detail::VecBase<T,3>::x+w*detail::VecBase<T,3>::z))*c+rsin*( w*detail::VecBase<T,3>::x-u*detail::VecBase<T,3>::z))/uvw2,
                    (w*uxvywz+(detail::VecBase<T,3>::z*(u2 + v2)-w*(u*detail::VecBase<T,3>::x+v*detail::VecBase<T,3>::y))*c+rsin*(-v*detail::VecBase<T,3>::x+u*detail::VecBase<T,3>::y))/uvw2);
        }
        
        const Vec<T,3> rotate(Vec<T,3> axis, Vec<T,3> center, T radians){
            return rotate(axis.x, axis.y, axis.z, center.x, center.y, center.z, radians);
        }
        
        const Vec<T,3> rotate(T u, T v, T w, T a, T b, T c, T radians){
            //rotation about arbitrary axis <u,v,w> and point <a,b,c>
            T x = detail::VecBase<T,3>::x;
            T y = detail::VecBase<T,3>::y;
            T z = detail::VecBase<T,3>::z;
            
            T cc = cos(radians);
            T ss = sin(radians);
            T oneMcos = 1 - cc;
            
            T u2 = u*u;
            T v2 = v*v;
            T w2 = w*w;
            
            T uxvywz = u*x - v*y - w*z;
            
            return Vec<T,3>(
                    (a*(v2 + w2) - u*(b*v + c*w - uxvywz))*oneMcos + x*cc + (-c*v + b*w - w*y + v*z)*ss,
                    (b*(u2 + w2) - v*(a*u + c*w - uxvywz))*oneMcos + y*cc + ( c*u - a*w + w*x - u*z)*ss,
                    (c*(u2 + v2) - w*(a*u + b*v - uxvywz))*oneMcos + z*cc + (-b*u + a*v - v*x + u*y)*ss);
        }
        
        #ifdef GEOMC_VEC_USE_SWIZZLE
        const Vec<T,3> swizzle(const std::string& xyz) const {
            if (xyz.length() != 3){
                throw GeomException("Invalid 3D swizzle string: " + xyz);
            }
            return Vec<T,3>(swizzlecoord(xyz[0]),
                            swizzlecoord(xyz[1]),
                            swizzlecoord(xyz[2]));
        }
        #endif

        bool isFiniteReal() const {
            return isReal(detail::VecBase<T,3>::x) && isReal(detail::VecBase<T,3>::y) && isReal(detail::VecBase<T,3>::z);
        }

    private:
        
        #ifdef GEOMC_VEC_USE_SWIZZLE
        T swizzlecoord(char c) const {
            switch (tolower(c)){
                case 'x':
                case 'r':
                    return detail::VecBase<T,3>::x;
                case 'y':
                case 'g':
                    return detail::VecBase<T,3>::y;
                case 'z':
                case 'b':
                    return detail::VecBase<T,3>::z;
                case '1':
                    return 1;
                case '0':
                    return 0;
                default:
                    throw GeomException("Invalid swizzle character: " + c);
            }
        }
        #endif
        
    }; //end Vec3 definition

    template <typename T> const Vec<T,3> Vec<T,3>::X_AXIS = Vec<T,3>(1,0,0);
    template <typename T> const Vec<T,3> Vec<T,3>::Y_AXIS = Vec<T,3>(0,1,0);
    template <typename T> const Vec<T,3> Vec<T,3>::Z_AXIS = Vec<T,3>(0,0,1);

} //end namespace geom



#endif /* VECTOR3D_H_ */
