/*
 * Vec3.h
 *
 *  Created on: Feb 22, 2009
 *      Author: Tim Babb
 */

#ifndef VECTOR3D_H_
#define VECTOR3D_H_

#include <cctype>
#include <geomc/linalg/vecdetail/VecBase.h>


namespace geom {

    /*===========================*
     * Global Functions          *
     *===========================*/

    template <typename T> 
    inline bool isReal(T d){
        return !(std::isinf(d) || std::isnan(d));
    }

    /**
     * @return The normal of the triangle formed by `p1`, `p2`, and `p3`.
     */
    template <typename T> 
    Vec<T,3> triangle_normal(const Vec<T,3> p1, const Vec<T,3> p2, const Vec<T,3> p3){
        return (p2-p1).cross(p3-p1).unit();
    }

    /**
     * An `(r,g,b)` color tuple representing a hue selected by `s`, with 0 being red,
     * progressing in rainbow order through the colors before returning to red at 1.
     * The cycle is extended similarly outside of `(0, 1)`.
     * 
     * @param s Rainbow parameter between 0 and 1. 
     * @return An `(r,g,b)` color tuple.
     */
    template <typename T> 
    Vec<T,3> rainbow(double s){
        s = std::fmod(s,1);
        if (s < 0) s = 1+s;
        return Vec<T,3>(
                std::min( std::max( std::max( -6*s+2, 6*s-4 ), 0.0), 1.0),
                std::max( std::min( std::min( 6*s,   -6*s+4 ), 1.0), 0.0),
                std::max( std::min( std::min( 6*s-2, -6*s+6 ), 1.0), 0.0));
    }

    /**
     * @param aRGB A packed 8-bit per channel RGB color value (with alpha 
     * or irrelevant data in the eight high order bits).
     * @return An `(r,g,b)` tuple extracted from `aRGB`.
     */
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

    /** @ingroup linalg
     *  @brief 3D specialization of vector class.
     * 
     * `Vec<T,3>`'s elements may be accessed under these equivalent naming schemes:
     * 
     *     v.{x,y,z} // conventional Euclidean coordinate names
     *     v.{s,t,u} // conventional parameterization coordinate names
     *     v.{r,g,b} // conventional color tuple coordinate names
     * 
     * `x`, `s`, and `r` all refer to the same element.
     */
    template <typename T> class Vec<T,3> : public detail::VecCommon< T, 3, Vec<T,3> > {
    public:
        
        /// @brief (1, 0, 0) constant
        static const Vec<T,3> X_AXIS;
        /// @brief (0, 1, 0) constant
        static const Vec<T,3> Y_AXIS;
        /// @brief (0, 0, 1) constant
        static const Vec<T,3> Z_AXIS;

        /*===========================*
         * Structors                 *
         *===========================*/

        /// Construct vector with all elements set to 0.
        Vec():detail::VecCommon< T, 3, Vec<T,3> >(){}
        
        /// Construct a vector with `(x, y)` taken from the 2D vector `xy`, and `z` as the final element.
        Vec(Vec<T,2> xy, T z){
            detail::VecBase<T,3>::x = xy.x;
            detail::VecBase<T,3>::y = xy.y;
            detail::VecBase<T,3>::z = z;
        }

        /// Construct a vector with all elements set to `a`.
        Vec(T a):detail::VecCommon< T, 3, Vec<T,3> >(a){}

        /// Construct a vector with elements copied from the 3-element array `v`.
        Vec(const T v[3]):detail::VecCommon< T, 3, Vec<T,3> >(v){}
        
        /// Construct a vector with elements `(x, y, z)`.
        Vec(T x, T y, T z){
            detail::VecBase<T,3>::x = x;
            detail::VecBase<T,3>::y = y;
            detail::VecBase<T,3>::z = z;
        }

        /*===========================*
         * Operators                 *
         *===========================*/

        /// @brief Vector cross product
        Vec<T,3>  operator^(Vec<T,3> v) const {
            return cross(v);
        }
        
        /*===========================*
         * Arithmetic                *
         *===========================*/
        
        //do not hide these; we want to overload them: 
        using detail::VecCommon< T, 3, Vec<T,3> >::add;
        using detail::VecCommon< T, 3, Vec<T,3> >::sub;
        using detail::VecCommon< T, 3, Vec<T,3> >::scale;
        using detail::VecCommon< T, 3, Vec<T,3> >::dot;
        
        /// @return `(dx, dy, dz)` added with `this`
        Vec<T,3> add(T dx, T dy, T dz) const {
            return Vec<T,3>(detail::VecBase<T,3>::x+dx, detail::VecBase<T,3>::y+dy, detail::VecBase<T,3>::z+dz);
        }

        /// @return `(dx, dy, dz)` subtracted from `this`
        Vec<T,3> sub(T dx, T dy, T dz) const {
            return Vec<T,3>(detail::VecBase<T,3>::x-dx, detail::VecBase<T,3>::y-dy, detail::VecBase<T,3>::z-dz);
        }

        /// @return Element-wise scaled copy
        Vec<T,3> scale(T a, T b, T c) const {
            return Vec<T,3>(a*detail::VecBase<T,3>::x, b*detail::VecBase<T,3>::y, c*detail::VecBase<T,3>::z);
        }
        
        /// @return Dot product of `this` with the vector `(xx, yy, zz)`
        T dot(T xx, T yy, T zz) const {
            return (detail::VecBase<T,3>::x * xx) + (detail::VecBase<T,3>::y*yy) + (detail::VecBase<T,3>::z*zz);
        }

        /// @return A copy with the X component negated
        Vec<T,3> reflectX() const {
            return Vec<T,3>(-detail::VecBase<T,3>::x, detail::VecBase<T,3>::y, detail::VecBase<T,3>::z);
        }

        /// @return A copy with the Y component negated
        Vec<T,3> reflectY() const {
            return Vec<T,3>(detail::VecBase<T,3>::x, -detail::VecBase<T,3>::y, detail::VecBase<T,3>::z);
        }

        /// @return A copy with the Z component negated
        Vec<T,3> reflectZ() const {
            return Vec<T,3>(detail::VecBase<T,3>::x, detail::VecBase<T,3>::y, -detail::VecBase<T,3>::z);
        }
        
        /*===========================*
         * 3D Geometry               *
         *===========================*/
        
        /// @return Cross product of `this` by `v`, in the order shown.
        Vec<T,3> cross(Vec<T,3> v) const {
            T newX = (detail::VecBase<T,3>::y * v.z) - (detail::VecBase<T,3>::z * v.y);
            T newY = (detail::VecBase<T,3>::z * v.x) - (detail::VecBase<T,3>::x * v.z);
            T newZ = (detail::VecBase<T,3>::x * v.y) - (detail::VecBase<T,3>::y * v.x);

            return Vec<T,3>(newX, newY, newZ);
        }
        
        /**
         * Convert this euclidean point to spherical coordinates `(r, elevation, azimuth)`,
         * with both angles in radians. Elevation is the angle from the +z axis in the 
         * range `(0, pi)`, and azimuth is the angle about +z from the +x axis
         * in the range `(-pi, pi)`.
         * 
         * @return This euclidean point converted to spherical coordinates.
         */
        Vec<T,3> toSpherical() const {
            T r = this->mag();
            return Vec<T,3>(r, std::acos(detail::VecBase<T,3>::z/r), std::atan2(detail::VecBase<T,3>::y,detail::VecBase<T,3>::x));
        }
        
        /**
         * Convert this point, interpreted as spherical coordinates `(r, elevation, azimuth)`,
         * to euclidean coordinates, where `r` (`this[0]`) is the radius, `eleveation`
         * (`this[1]`) is the angle from the +z axis, and `azimuth` (`this[2]`)
         * is the angle about z+ from the x+ axis.
         * 
         * @return This point described by this vector interpreted as spherical 
         * coordinates as `(x, y, z)` euclidean coordinates.
         */
        Vec<T,3> fromSpherical() const {
            Vec<T,3> oot;
            T rsinelev = detail::VecBase<T,3>::x * std::sin(detail::VecBase<T,3>::y);
            oot.x = rsinelev*std::cos(detail::VecBase<T,3>::z);
            oot.y = rsinelev*std::sin(detail::VecBase<T,3>::z);
            oot.z = detail::VecBase<T,3>::x*std::cos(detail::VecBase<T,3>::y);
            return oot;
        }

        /// @return A copy of this vector rotated by angle `radians` about `axis`.
        Vec<T,3> rotate(Vec<T,3> axis, T radians) const {
            return rotate(axis.x, axis.y, axis.z, radians);
        }

        /// @return A copy of this vector rotated by angle `radians` about the axis `(u, v, w)`.
        Vec<T,3> rotate(T u, T v, T w, T radians) const {
            T u2 = u*u;
            T v2 = v*v;
            T w2 = w*w;
            T uxvywz = u*detail::VecBase<T,3>::x+v*detail::VecBase<T,3>::y+w*detail::VecBase<T,3>::z;
            T uvw2 = u2+v2+w2;
            T rsin = std::sqrt(uvw2)*std::sin(radians);
            T c    = std::cos(radians);

            return Vec<T,3>(
                    (u*uxvywz+(detail::VecBase<T,3>::x*(v2 + w2)-u*(v*detail::VecBase<T,3>::y+w*detail::VecBase<T,3>::z))*c+rsin*(-w*detail::VecBase<T,3>::y+v*detail::VecBase<T,3>::z))/uvw2,
                    (v*uxvywz+(detail::VecBase<T,3>::y*(u2 + w2)-v*(u*detail::VecBase<T,3>::x+w*detail::VecBase<T,3>::z))*c+rsin*( w*detail::VecBase<T,3>::x-u*detail::VecBase<T,3>::z))/uvw2,
                    (w*uxvywz+(detail::VecBase<T,3>::z*(u2 + v2)-w*(u*detail::VecBase<T,3>::x+v*detail::VecBase<T,3>::y))*c+rsin*(-v*detail::VecBase<T,3>::x+u*detail::VecBase<T,3>::y))/uvw2);
        }
        
        /** @return A copy of this vector rotated by angle `radians` about 
         * rotation axis `axis` and centerpoint `center`.
         */
        Vec<T,3> rotate(Vec<T,3> axis, Vec<T,3> center, T radians){
            return rotate(axis.x, axis.y, axis.z, center.x, center.y, center.z, radians);
        }
        
        /**
         * @return A copy of this vector rotated by angle `radians` about
         * the axis `(u, v, w)`, and centerpoint `(a, b, c)`.
         */
        Vec<T,3> rotate(T u, T v, T w, T a, T b, T c, T radians){
            T x = detail::VecBase<T,3>::x;
            T y = detail::VecBase<T,3>::y;
            T z = detail::VecBase<T,3>::z;
            
            T cc = std::cos(radians);
            T ss = std::sin(radians);
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

        /// @return `true` if no elements are an infinity or `NaN`.
        bool isFiniteReal() const {
            return isReal(detail::VecBase<T,3>::x) && isReal(detail::VecBase<T,3>::y) && isReal(detail::VecBase<T,3>::z);
        }

    private:
        
    }; //end Vec3 definition

    template <typename T> const Vec<T,3> Vec<T,3>::X_AXIS = Vec<T,3>(1,0,0);
    template <typename T> const Vec<T,3> Vec<T,3>::Y_AXIS = Vec<T,3>(0,1,0);
    template <typename T> const Vec<T,3> Vec<T,3>::Z_AXIS = Vec<T,3>(0,0,1);

} //end namespace geom



#endif /* VECTOR3D_H_ */
