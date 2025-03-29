#pragma once

/*
 * Vec3.h
 *
 *  Created on: Feb 22, 2009
 *      Author: Tim Babb
 */

#include <cctype>
#include <geomc/linalg/vecdetail/VecBase.h>
#include <geomc/function/Utils.h>


namespace geom {

    /*===========================*
     * Global Functions          *
     *===========================*/

    /**
     * @return The normal of the triangle formed by `p1`, `p2`, and `p3`.
     */
    template <typename T> 
    Vec<T,3> triangle_normal(const Vec<T,3> p1, const Vec<T,3> p2, const Vec<T,3> p3){
        return (p2 - p1).cross(p3 - p1).unit();
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
    Vec<T,3> from_rgb(int aRGB){
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
    
        using base_t = detail::VecCommon<T,3,Vec<T,3>>;

        // to prevent c++20 ambiguity warnings:
        using detail::VecCommon<T,3,Vec<T,3>>::operator==;
        
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
        constexpr Vec() {}
        
        /**
         * @brief Construct a vector with `(x, y)` taken from the 2D vector `xy`,
         * and `z` as the final element.
         */
        template <typename U>
        constexpr Vec(Vec<U,2> xy, T z) {
            detail::VecBase<T,3>::x = xy.x;
            detail::VecBase<T,3>::y = xy.y;
            detail::VecBase<T,3>::z = z;
        }

        /// Construct a vector with all elements set to `a`.
        constexpr Vec(T a):base_t(a) {}

        /// Construct a vector with elements copied from the 3-element array `v`.
        constexpr Vec(const T v[3]):base_t(v) {}
        
        /// Construct a vector with elements `(x, y, z)`.
        constexpr Vec(T x, T y, T z) {
            detail::VecBase<T,3>::x = x;
            detail::VecBase<T,3>::y = y;
            detail::VecBase<T,3>::z = z;
        }

        /*===========================*
         * Operators                 *
         *===========================*/

        /// @brief Vector cross product
        Vec<T,3> operator^(Vec<T,3> v) const {
            return cross(v);
        }
        
        /*===========================*
         * Arithmetic                *
         *===========================*/
        
        // do not hide these; we want to overload them: 
        using base_t::add;
        using base_t::sub;
        using base_t::scale;
        using base_t::dot;
        
        /// @return `(dx, dy, dz)` added with `this`
        Vec<T,3> add(T dx, T dy, T dz) const {
            return Vec<T,3>(
                this->x + dx,
                this->y + dy,
                this->z + dz
            );
        }

        /// @return `(dx, dy, dz)` subtracted from `this`
        Vec<T,3> sub(T dx, T dy, T dz) const {
            return Vec<T,3>(
                this->x - dx,
                this->y - dy,
                this->z - dz
            );
        }

        /// @return Element-wise scaled copy
        Vec<T,3> scale(T a, T b, T c) const {
            return Vec<T,3>(
                a*this->x,
                b*this->y,
                c*this->z
            );
        }
        
        /// @return Dot product of `this` with the vector `(xx, yy, zz)`
        T dot(T xx, T yy, T zz) const {
            return
                (this->x * xx) +
                (this->y * yy) +
                (this->z * zz);
        }

        /// @return A copy with the X component negated
        inline Vec<T,3> reflected_x() const {
            return Vec<T,3>(
                -this->x,
                 this->y,
                 this->z
            );
        }

        /// @return A copy with the Y component negated
        inline Vec<T,3> reflected_y() const {
            return Vec<T,3>(
                 this->x,
                -this->y,
                 this->z
            );
        }

        /// @return A copy with the Z component negated
        inline Vec<T,3> reflected_z() const {
            return Vec<T,3>(
                 this->x,
                 this->y,
                -this->z
            );
        }
        
        /*===========================*
         * 3D Geometry               *
         *===========================*/
        
        /// @return Cross product of `this` by `v`, in the order shown.
        Vec<T,3> cross(Vec<T,3> v) const {
            T newX = diff_of_products(
                this->y, v.z,
                this->z, v.y);
            T newY = diff_of_products(
                this->z, v.x,
                this->x, v.z);
            T newZ = diff_of_products(
                this->x, v.y,
                this->y, v.x);

            return Vec<T,3>(newX, newY, newZ);
        }
        
        /**
         * Convert this euclidean point to spherical coordinates `(r, polar, azimuth)`,
         * with both angles in radians. Polar angle is the angle from the +z axis in the 
         * range `(0, pi)`, and azimuth is the angle about +z from the +x axis
         * in the range `(-pi, pi)`.
         * 
         * @return This euclidean point converted to spherical coordinates.
         */
        Vec<T,3> to_spherical() const {
            T r = this->mag();
            return Vec<T,3>(
                r,
                std::acos(detail::VecBase<T,3>::z/r),
                std::atan2(
                    detail::VecBase<T,3>::y,
                    detail::VecBase<T,3>::x));
        }
        
        /**
         * Construct a point from spherical coordinates `(r, polar, azimuth)`,
         * to euclidean coordinates, where `r` is the radius, `polar`
         * is the angle from the +z axis, and `azimuth` is the angle about z+ from
         * the x+ axis.
         * 
         * @return A euclidean point constructed from spherical coordinates.
         */
        static Vec<T,3> from_spherical(T r, T polar, T azi) {
            Vec<T,3> out;
            T rsinpolar = r * std::sin(polar);
            out.x = rsinpolar * std::cos(azi);
            out.y = rsinpolar * std::sin(azi);
            out.z = r * std::cos(polar);
            return out;
        }

        /// @return A copy of this vector rotated by angle `radians` about `axis`.
        Vec<T,3> rotate(Vec<T,3> axis, T radians) const {
            return rotate(axis.x, axis.y, axis.z, radians);
        }

        /**
         * @return A copy of this vector rotated by angle `radians` about
         * the axis `(u, v, w)`.
         */
        Vec<T,3> rotate(T u, T v, T w, T radians) const {
            const T& x = this->x;
            const T& y = this->y;
            const T& z = this->z;
            T u2 = u * u;
            T v2 = v * v;
            T w2 = w * w;  // likely to become an FMA, which will keep precision ↙︎
            T uxvywz = sum_of_products(u, x, v, y) + w * z;
            T uvw2 = u2 + v2 + w2;
            T rsin = std::sqrt(uvw2) * std::sin(radians);
            T c    = std::cos(radians);
            
            // intermediate sums, using a safe high-precision formula
            // (expected to be within 5% of "naïve" performance)
            T p0 = sum_of_products(v, y, w, z);
            T p1 = sum_of_products(u, x, w, z);
            T p2 = sum_of_products(u, z, v, y);
            
            T q0 = diff_of_products(v, z, w, y);
            T q1 = diff_of_products(w, x, u, z);
            T q2 = diff_of_products(u, y, v, x);
            
            T s0 = diff_of_products(x, v2 + w2, u, p0);
            T s1 = diff_of_products(y, u2 + w2, v, p1);
            T s2 = diff_of_products(z, u2 + v2, w, p2);
            
            T t0 = sum_of_products(u, uxvywz, s0, c);
            T t1 = sum_of_products(v, uxvywz, s1, c);
            T t2 = sum_of_products(w, uxvywz, s2, c);
            
            return Vec<T,3>(
                    (t0 + rsin * q0) / uvw2,  // <- also likely to become FMA
                    (t1 + rsin * q1) / uvw2,
                    (t2 + rsin * q2) / uvw2);
        }
        
        /** @return A copy of this vector rotated by angle `radians` about 
         * rotation axis `axis` and centerpoint `center`.
         */
        inline Vec<T,3> rotate(Vec<T,3> axis, Vec<T,3> center, T radians){
            return rotate(axis.x, axis.y, axis.z, center.x, center.y, center.z, radians);
        }
        
        /**
         * @return A copy of this vector rotated by angle `radians` about
         * the axis `(u, v, w)`, and centerpoint `(a, b, c)`.
         */
        Vec<T,3> rotate(T u, T v, T w, T a, T b, T c, T radians){
            const T& x = this->x;
            const T& y = this->y;
            const T& z = this->z;
            
            T cc = std::cos(radians);
            T ss = std::sin(radians);
            T oneMcos = 1 - cc;
            
            T u2 = u * u;
            T v2 = v * v;
            T w2 = w * w;
            T uxvywz = diff_of_products(u, x, v, y) - (w * z);
            
            T p0 = sum_of_products(b, v, c, w);
            T p1 = sum_of_products(a, u, c, w);
            T p2 = sum_of_products(a, u, b, v);
            T q0 = diff_of_products(a, v2 + w2, u, p0 - uxvywz);
            T q1 = diff_of_products(b, u2 + w2, v, p1 - uxvywz);
            T q2 = diff_of_products(c, u2 + v2, w, p2 - uxvywz);
            
            return Vec<T,3>(
                    q0 * oneMcos + x * cc + (-c*v + b*w - w*y + v*z) * ss,
                    q1 * oneMcos + y * cc + ( c*u - a*w + w*x - u*z) * ss,
                    q2 * oneMcos + z * cc + (-b*u + a*v - v*x + u*y) * ss);
        }
        
        /**
         * Compute a transform aligning `from` with `to` and apply it to `this`.
         */
        Vec<T,3> align(const Vec<T,3>& from, const Vec<T,3>& to) {
            // modified from http://iquilezles.org/www/articles/noacos/noacos.htm
            const Vec<T,3>& v = *this;
            const Vec<T,3>  q = from.cross(to);
            const T         c = from.dot(to);
            const T         k = 1 / (1 + c);
            return Vec<T,3>(
                (q.x * q.x * k + c)   * v.x + (q.y * q.x * k - q.z) * v.y + (q.z * q.x * k + q.y) * v.z,
                (q.x * q.y * k + q.z) * v.x + (q.y * q.y * k + c)   * v.y + (q.z * q.y * k - q.x) * v.z,
                (q.x * q.z * k - q.y) * v.x + (q.y * q.z * k + q.x) * v.y + (q.z * q.z * k + c)   * v.z);
        }
        
        /// @return `true` if no elements are an infinity or `NaN`.
        bool is_finite_real() const {
            return is_real(this->x) and
                   is_real(this->y) and
                   is_real(this->z);
        }

    private:
        
    }; //end Vec3 definition

    template <typename T> const Vec<T,3> Vec<T,3>::X_AXIS = Vec<T,3>(1,0,0);
    template <typename T> const Vec<T,3> Vec<T,3>::Y_AXIS = Vec<T,3>(0,1,0);
    template <typename T> const Vec<T,3> Vec<T,3>::Z_AXIS = Vec<T,3>(0,0,1);

} //end namespace geom
