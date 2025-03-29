#pragma once
/*
 * Vec4.h
 *
 *  Created on: Feb 22, 2009
 *      Author: Tim Babb
 */

#include <cctype>
#include <algorithm>
#include <geomc/linalg/vecdetail/VecBase.h>

namespace geom {

    /*===========================*
     * Global Functions          *
     *===========================*/

    template <typename T> const Vec<T,4> from_argb(int aRGB) {
        T a = ((aRGB & 0xff000000) >> 24) / detail::_RGBChannelConversion<T>::factor;
        T r = ((aRGB & 0x00ff0000) >> 16) / detail::_RGBChannelConversion<T>::factor;
        T g = ((aRGB & 0x0000ff00) >> 8)  / detail::_RGBChannelConversion<T>::factor;
        T b =  (aRGB & 0x000000ff)        / detail::_RGBChannelConversion<T>::factor;
        return Vec<T,4>(r,g,b,a);
    }

    template <typename T> const Vec<T,4> from_rgba(int rgba) {
        T r = ((rgba & 0xff000000) >> 24) / detail::_RGBChannelConversion<T>::factor;
        T g = ((rgba & 0x00ff0000) >> 16) / detail::_RGBChannelConversion<T>::factor;
        T b = ((rgba & 0x0000ff00) >> 8)  / detail::_RGBChannelConversion<T>::factor;
        T a =  (rgba & 0x000000ff)        / detail::_RGBChannelConversion<T>::factor;
        return Vec<T,4>(r,g,b,a);
    }
    
    /**
     * @brief Composite two non-premultiplied colors over each other.
     */
    template <typename T>
    const Vec<T,4> composite_over(Vec<T,4> fg, Vec<T,4> bg) {
        Vec<T,3> c_fg = fg.template resized<3>();
        Vec<T,3> c_bg = bg.template resized<3>();
        T a_fg = fg.w;
        T a_bg = bg.w;
        T a_out = a_fg + (1 - a_fg) * a_bg;
        if (a_out == 0) {
            return fg + bg;
        }
        return {
            (c_fg * a_fg + c_bg * a_bg * (1 - a_fg)) / a_out,
            a_out
        };
    }

    /*===========================*
     * Vector Class              *
     *===========================*/

    /** @ingroup linalg
     *  @brief 4D specialization of vector class.
     * 
     * `Vec<T,4>`'s elements may be accessed under these equivalent naming schemes:
     * 
     *     v.{x,y,z,w} // conventional Euclidean coordinate names
     *     v.{s,t,u,v} // conventional parameterization coordinate names
     *     v.{r,g,b,a} // conventional color tuple coordinate names
     * 
     * with `x`, `s`, and `r` all referring to the same element.
     */
    template <typename T> class Vec<T,4> : public detail::VecCommon< T, 4, Vec<T,4> > {
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
        constexpr Vec():detail::VecCommon<T, 4, Vec<T,4>>() {}
        
        /// Construct a vector with elements `(x, y, z, w)`.
        constexpr Vec(T x, T y, T z, T w):detail::VecCommon<T, 4, Vec<T,4>>{x,y,z,w} {}

        /// Construct a vector with elements copied from the 4-element array `v`.
        constexpr Vec(const T v[4]):detail::VecCommon<T, 4, Vec<T,4>>(v) {};

        /**
         * @brief Construct a vector with `(x, y, z)` taken from the 3D vector `xyz`,
         * and `w` as the final element.
         */
        template <typename U>
        constexpr Vec(Vec<U,3> xyz, T w):detail::VecCommon<T, 4, Vec<T,4>>{
            xyz.x, xyz.y, xyz.z, w
        } {}
        
        template <typename U>
        constexpr Vec(Vec<U,2> xy, T z, T w):detail::VecCommon<T, 4, Vec<T,4>>{
            xy.x, xy.y, z, w
        } {}

        /// Construct a vector with all elements set to `a`.
        constexpr Vec(T a):detail::VecCommon<T, 4, Vec<T,4>>(a) {}
        
        /*===========================*
         * Arithmetic                *
         *===========================*/

        // to prevent c++20 ambiguity warnings:
        using detail::VecCommon<T,4,Vec<T,4>>::operator==;
        
        using detail::VecCommon< T, 4, Vec<T,4> >::add;
        using detail::VecCommon< T, 4, Vec<T,4> >::sub;
        using detail::VecCommon< T, 4, Vec<T,4> >::scale;
        using detail::VecCommon< T, 4, Vec<T,4> >::dot;

        /// @return `(dx, dy, dz, dw)` added with `this`
        Vec<T,4> add(T dx, T dy, T dz, T dw) const {
            return Vec<T,4>(detail::VecBase<T,4>::x+dx, 
                            detail::VecBase<T,4>::y+dy, 
                            detail::VecBase<T,4>::z+dz, 
                            detail::VecBase<T,4>::w+dw);
        }
        
        /// @return `(dx, dy, dz, dw)` subtracted from `this`
        Vec<T,4> sub(T dx, T dy, T dz, T dw) const {
            return Vec<T,4>(detail::VecBase<T,4>::x-dx, 
                            detail::VecBase<T,4>::y-dy, 
                            detail::VecBase<T,4>::z-dz, 
                            detail::VecBase<T,4>::w-dw);
        }
        
        /// @return Element-wise scaled copy
        Vec<T,4> scale(T a, T b, T c, T d) const {
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
        
    }; //end Vec4 definition
    
    template <typename T> const Vec<T,4> Vec<T,4>::X_AXIS = Vec<T,4>(1,0,0,0);
    template <typename T> const Vec<T,4> Vec<T,4>::Y_AXIS = Vec<T,4>(0,1,0,0);
    template <typename T> const Vec<T,4> Vec<T,4>::Z_AXIS = Vec<T,4>(0,0,1,0);
    template <typename T> const Vec<T,4> Vec<T,4>::W_AXIS = Vec<T,4>(0,0,0,1);

} //end namespace geom
