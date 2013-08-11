/*
 * Vec4.h
 *
 *  Created on: Feb 22, 2009
 *      Author: Tim Babb
 */

#ifndef VECTOR4D_H_
#define VECTOR4D_H_

#include <cctype>
#include <algorithm>
#include <geomc/linalg/vecdetail/VecBase.h>

#ifdef GEOMC_VEC_USE_SWIZZLE
    #include <string>
    #include <geomc/GeomException.h>
#endif

namespace geom {

    /*===========================*
     * Global Functions          *
     *===========================*/

    template <typename T> const Vec<T,4> fromARGB(int aRGB){
        T a = ((aRGB & 0xff000000) >> 24) / detail::_RGBChannelConversion<T>::factor;
        T r = ((aRGB & 0x00ff0000) >> 16) / detail::_RGBChannelConversion<T>::factor;
        T g = ((aRGB & 0x0000ff00) >> 8)  / detail::_RGBChannelConversion<T>::factor;
        T b =  (aRGB & 0x000000ff)        / detail::_RGBChannelConversion<T>::factor;
        return Vec<T,4>(r,g,b,a);
    }

    template <typename T> const Vec<T,4> fromRGBA(int rgba){
        T r = ((rgba & 0xff000000) >> 24) / detail::_RGBChannelConversion<T>::factor;
        T g = ((rgba & 0x00ff0000) >> 16) / detail::_RGBChannelConversion<T>::factor;
        T b = ((rgba & 0x0000ff00) >> 8)  / detail::_RGBChannelConversion<T>::factor;
        T a =  (rgba & 0x000000ff)        / detail::_RGBChannelConversion<T>::factor;
        return Vec<T,4>(r,g,b,a);
    }

    /*===========================*
     * Vector Class              *
     *===========================*/

    /** @ingroup linalg
     *  @brief 4D specialization of vector class.
     * 
     * `Vec<T,4>`'s elements may be accessed under these equivalent naming schemes:
     * 
     *     v.{x,y,z,w}
     *     v.{s,t,u,v}
     *     v.{r,g,b,a}
     * 
     * with `x`, `s`, and `r` all referring to the same element.
     */
    template <typename T> class Vec<T,4> : public detail::VecCommon<T,4> {
    public:
        
        /// @brief (1, 0, 0, 0) constant
        static const Vec<T,4> X_AXIS;
        /// @brief (0, 1, 0, 0) constant
        static const Vec<T,4> Y_AXIS;
        /// @brief (0, 0, 1, 0) constant
        static const Vec<T,4> Z_AXIS;
        /// @brief (0, 0, 0, 1) constant
        static const Vec<T,4> W_AXIS;

        /*===========================*
         * Structors                 *
         *===========================*/

        
        /// Construct vector with all elements set to 0.
        Vec():detail::VecCommon<T,4>(){}
        
        /// Construct a vector with elements `(x, y, z, w)`.
        Vec(T x, T y, T z, T w){
            detail::VecBase<T,4>::x = x;
            detail::VecBase<T,4>::y = y;
            detail::VecBase<T,4>::z = z;
            detail::VecBase<T,4>::w = w;
        }

        /// Construct a vector with elements copied from the 4-element array `v`.
        Vec(T v[4]):detail::VecCommon<T,4>(v){};

        /// Construct a vector with `(x, y, z)` taken from the 3D vector `xyz`, and `w` as the final element.
        Vec(Vec<T,3> xyz, T w){
            detail::VecBase<T,4>::x = xyz.x;
            detail::VecBase<T,4>::y = xyz.y;
            detail::VecBase<T,4>::z = xyz.z;
            detail::VecBase<T,4>::w = w;
        }

        /// Construct a vector with all elements set to `a`.
        Vec(T a):detail::VecCommon<T,4>(a){}
        
        /*===========================*
         * Arithmetic                *
         *===========================*/
        
        using detail::VecCommon<T,4>::add;
        using detail::VecCommon<T,4>::sub;
        using detail::VecCommon<T,4>::scale;
        using detail::VecCommon<T,4>::dot;

        /// @return `(dx, dy, dz, dw)` added with `this`
        const Vec<T,4> add(T dx, T dy, T dz, T dw) const {
            return Vec<T,4>(detail::VecBase<T,4>::x+dx, 
                            detail::VecBase<T,4>::y+dy, 
                            detail::VecBase<T,4>::z+dz, 
                            detail::VecBase<T,4>::w+dw);
        }
        
        /// @return `(dx, dy, dz, dw)` subtracted from `this`
        const Vec<T,4> sub(T dx, T dy, T dz, T dw) const {
            return Vec<T,4>(detail::VecBase<T,4>::x-dx, 
                            detail::VecBase<T,4>::y-dy, 
                            detail::VecBase<T,4>::z-dz, 
                            detail::VecBase<T,4>::w-dw);
        }
        
        /// @return Element-wise scaled copy
        const Vec<T,4> scale(T a, T b, T c, T d) const {
            return Vec<T,4>(a*detail::VecBase<T,4>::x, 
                            b*detail::VecBase<T,4>::y, 
                            c*detail::VecBase<T,4>::z, 
                            d*detail::VecBase<T,4>::w);
        }

        /// @return Dot product of `this` with the vector `(xx, yy, zz, ww)`
        T dot(T xx, T yy, T zz, T ww) const {
            return (detail::VecBase<T,4>::x * xx) + 
                   (detail::VecBase<T,4>::y * yy) + 
                   (detail::VecBase<T,4>::z * zz) + 
                   (detail::VecBase<T,4>::w * ww);
        }
        
        /*===========================*
         * 4D Geometry               *
         *===========================*/

        #ifdef GEOMC_VEC_USE_SWIZZLE
        const Vec<T,4> swizzle(const std::string& xyzw) const {
            if (xyzw.length() != 4){
                throw GeomException("Invalid 4D swizzle string: " + xyzw);
            }
            return Vec<T,4>(swizzlecoord(xyzw[0]),
                              swizzlecoord(xyzw[1]),
                            swizzlecoord(xyzw[2]),
                            swizzlecoord(xyzw[3]));
        }
        #endif

    private:
        
        #ifdef GEOMC_VEC_USE_SWIZZLE
        inline T swizzlecoord(char c) const {
            switch (tolower(c)){
                case 'x':
                case 'r':
                    return detail::VecBase<T,4>::x;
                case 'y':
                case 'g':
                    return detail::VecBase<T,4>::y;
                case 'z':
                case 'b':
                    return detail::VecBase<T,4>::z;
                case 'w':
                case 'a':
                    return detail::VecBase<T,4>::w;
                case '1':
                    return 1;
                case '0':
                    return 0;
                default:
                    throw GeomException("Invalid swizzle character: " + c);
            }
        }
        #endif
        
    }; //end Vec4 definition
    
    template <typename T> const Vec<T,4> Vec<T,4>::X_AXIS = Vec<T,4>(1,0,0,0);
    template <typename T> const Vec<T,4> Vec<T,4>::Y_AXIS = Vec<T,4>(0,1,0,0);
    template <typename T> const Vec<T,4> Vec<T,4>::Z_AXIS = Vec<T,4>(0,0,1,0);
    template <typename T> const Vec<T,4> Vec<T,4>::W_AXIS = Vec<T,4>(0,0,0,1);

} //end namespace geom

#endif /* VECTOR4D_H_ */
