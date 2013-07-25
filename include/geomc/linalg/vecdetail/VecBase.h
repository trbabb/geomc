/*
 * VecBase.h
 *
 *  Created on: Oct 10, 2010
 *      Author: tbabb
 */

#ifndef VECBASE_H_
#define VECBASE_H_

#include <cmath>
#include <algorithm>
#include <tr1/functional>
#include <boost/integer.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_integral.hpp>

#include <geomc/linalg/LinalgTypes.h>

#ifdef GEOMC_VEC_CHECK_BOUNDS
    #include <stdexcept>
#endif

#ifdef GEOMC_LINALG_USE_STREAMS
    #include <iostream>
#endif

//TODO: specialize for a special dimension Vec<t,DYNAMIC=-1>

//TODO: use c++11's decltype() and <auto> functionality to do proper widening/narrowing conversion
//      automatically, per http://stackoverflow.com/questions/9368432/c-arithmetic-operator-overloadingautomatic-widening
//      see also: boost::common_type<T,U>

//TODO: use begin() and end() pointers and pointer iteration instead of a counter?

namespace geom {
    
namespace detail {

    template <typename T, index_t N> class VecBase {
    public:

        VecBase(){
            std::fill(begin(), end(), 0);
        }

        VecBase(T a){
            std::fill(begin(), end(), a);
        }

        VecBase(T a[N]){
            std::copy(a, a+N, begin());
        }
        
        inline const T& get(index_t idx) const {
            #ifdef GEOMC_VEC_CHECK_BOUNDS
                // this is considerably slower
                if (idx < 0 or idx >= N){
                    throw std::out_of_range("vector index");
                }
            #endif
            return v[idx];
        }
        
        inline T& get(index_t idx){
            #ifdef GEOMC_VEC_CHECK_BOUNDS
                if (idx < 0 or idx >= N){
                    throw std::out_of_range("vector index");
                }
            #endif
            return v[idx];
        }
        
        inline const T* begin() const {
            return v;
        }
        
        inline const T* end() const {
            return v+N;
        }
        
        inline T* begin() {
            return v;
        }
        
        inline T* end() {
            return v+N;
        }
        
    protected:
        T v[N];
    };
    
    template <typename T> class VecBase<T,2> {
    public:

        VecBase():x(0),y(0){}

        VecBase(T a):x(a),y(a){}

        VecBase(T a[2]):x(a[0]),y(a[1]){}
        
        inline const T& get(index_t idx) const {
            #ifdef GEOMC_VEC_CHECK_BOUNDS
            switch(idx){
            case 0:
                return x;
            case 1:
                return y;
            default:
                throw std::out_of_range("vector index");
            }
            #else
            return idx ? y : x;
            #endif
        }
        
        inline T& get(index_t idx) {
            #ifdef GEOMC_VEC_CHECK_BOUNDS
            switch(idx){
            case 0:
                return x;
            case 1:
                return y;
            default:
                throw std::out_of_range("vector index");
            }
            #else
                return idx ? y : x;
            #endif
        }
        
        inline const T* begin() const {
            return &x;
        }
        
        inline const T* end() const {
            return &y + 1;
        }
        
        inline T* begin() {
            return &x;
        }
        
        inline T* end() {
            return &y + 1;
        }
        
        union {
            T x;
            T s;
            T row; // for matrix coordinates
        };
        
        union {
            T y;
            T t;
            T col;
        };
    };
    
    template <typename T> class VecBase<T,3> {
    public:

        VecBase():x(0),y(0),z(0){}

        VecBase(T a):x(a),y(a),z(a){}

        VecBase(T a[3]):x(a[0]),y(a[1]),z(a[2]){}

        inline const T& get(index_t idx) const {
        #ifdef GEOMC_VEC_CHECK_BOUNDS
            switch (idx){
                case 0:
                    return x;
                case 1:
                    return y;
                case 2:
                    return z;
                default:
                    throw std::out_of_range("vector index"); //throw is very expensive, even if we don't use it. consider avoiding this.
            }
        #else
            return (&x)[idx]; //faster than a 4-branch case, whether the case has a throw or not (gcc)
        #endif
        }
        
        inline T& get(index_t idx) {
#ifdef GEOMC_CHECK_BOUNDS
            switch (idx){
                case 0:
                    return x;
                case 1:
                    return y;
                case 2:
                    return z;
                default:
                    throw std::out_of_range("vector index"); //very slow; whether used or not.
            }
#else
            return (&x)[idx];
#endif
        }
        
        inline const T* begin() const {
            return &x;
        }
        
        inline const T* end() const {
            return &z + 1;
        }
        
        inline T* begin() {
            return &x;
        }
        
        inline T* end() {
            return &z + 1;
        }
        
        union {
            T x;
            T s;
            T r;
        };
        union {
            T y;
            T t;
            T g;
        };
        union {
            T z;
            T u;
            T b;
        };
    };
    
    template <typename T> class VecBase<T,4> {
    public:

        VecBase():x(0),y(0),z(0),w(0){}

        VecBase(T a):x(a),y(a),z(a),w(a){}

        VecBase(T a[4]):x(a[0]),y(a[1]),z(a[2]),w(a[3]){}

        inline const T& get(index_t idx) const {
#ifdef GEOMC_VEC_CHECK_BOUNDS
            switch (idx){
                case 0:
                    return x;
                case 1:
                    return y;
                case 2:
                    return z;
                case 3:
                    return w;
                default:
                    throw std::out_of_range("vector index"); //very slow; whether used or not.
            }
#else
            return (&x)[idx];
#endif
        }
        
        inline T& get(index_t idx) {
#ifdef GEOMC_VEC_CHECK_BOUNDS
            switch (idx){
                case 0:
                    return x;
                case 1:
                    return y;
                case 2:
                    return z;
                case 3:
                    return w;
                default:
                    throw std::out_of_range("vector index");
            }
#else
            return (&x)[idx];
#endif
        }
        
        inline const T* begin() const {
            return &x;
        }
        
        inline const T* end() const {
            return &w + 1;
        }
        
        inline T* begin() {
            return &x;
        }
        
        inline T* end() {
            return &w + 1;
        }
        
        union {
            T x;
            T s;
            T r;
        };
        union {
            T y;
            T t;
            T g;
        };
        union {
            T z;
            T u;
            T b;
        };
        union {
            T w;
            T v;
            T a;
        };
    };
    
    // end VecBase classes

    template <typename T, index_t N> class VecCommon : public VecBase<T,N> {
        public:
        
        typedef T elem_t;
        static const index_t DIM = N;
        
        static const Vec<T,N> ones;
        static const Vec<T,N> zeros;

        /*******************************
         * Structors                   *
         *******************************/

        VecCommon():VecBase<T,N>(){}

        VecCommon(T a):VecBase<T,N>(a){}

        VecCommon(T a[N]):VecBase<T,N>(a){}

        /*******************************
         * Operators                   *
         *******************************/

        inline const T& operator[](index_t idx) const {
            return this->get(idx);
        }
        
        inline T& operator[](index_t idx) {
            return this->get(idx);
        }

        template <typename U> operator Vec<U,N>() const {
            Vec<U,N> r;
            for (index_t i = 0; i < N; i++){
                r[i] = (U)(this->get(i));
            }
            return r;
        }

        inline bool operator==(const Vec<T,N> &vv) const {
            for (index_t i = 0; i < N; i++){
                if (vv.get(i) != this->get(i)) return false;
            }
            return true;
        }
        
        inline bool operator!=(const Vec<T,N> &vv) const {
            for (index_t i = 0; i < N; i++){
                if (vv.get(i) != this->get(i)) return true;
            }
            return false;
        }

        inline const Vec<T,N> operator+(const Vec<T,N> &v) const {
            return this->add(v);
        }

        VecCommon<T,N>& operator+=(const Vec<T,N> &vv){
            for (index_t i = 0; i < N; i++){
                this->get(i) += vv[i];
            }
            return *this;
        }
        
        inline const Vec<T,N>  operator-(const Vec<T,N> &v) const {
            return this->minus(v);
        }

        VecCommon<T,N>& operator-=(const Vec<T,N> &vv){
            for (index_t i = 0; i < N; i++){
                this->get(i) -= vv[i];
            }
            return *this;
        }

        VecCommon<T,N>& operator*=(T s){
            for (index_t i = 0; i < N; i++){
                this->get(i) *= s;
            }
            return *this;
        }

        VecCommon<T,N>& operator/=(T s){
            for (index_t i = 0; i < N; i++){
                this->get(i) /= s;
            }
            return *this;
        }

        VecCommon<T,N>& operator*=(const Vec<T,N> &vv){
            for (index_t i = 0; i < N; i++){
                this->get(i) *= vv[i];
            }
            return *this;
        }

        //negation
        inline const Vec<T,N> operator-() const {
            return this->neg();
        }

        /*******************************
         * Methods                     *
         *******************************/
        
        inline const Vec<T,N> add(const Vec<T,N> &v) const {
            Vec<T,N> r;
            for (index_t i = 0; i < N; i++){
                r[i] = this->get(i) + v.get(i);
            }
            return r;
        }
        
        inline const Vec<T,N> minus(const Vec<T,N> &v) const {
            Vec<T,N> r;
            for (index_t i = 0; i < N; i++){
                r[i] = this->get(i) - v[i];
            }
            return r;
        }
        
        inline const Vec<T,N> scale(const Vec<T,N> &v) const {
            Vec<T,N> r;
            for (index_t i = 0; i < N; i++){
                r[i] = this->get(i) * v[i];
            }
            return r;
        }
        
        inline const Vec<T,N> scale(T a) const {
            Vec<T,N> r;
            for (index_t i = 0; i < N; i++){
                r[i] = this->get(i) * a;
            }
            return r;
        }
        
        inline const Vec<T,N> neg() const {
            Vec<T,N> r;
            for (index_t i = 0; i < N; i++){
                r[i] = -this->get(i);
            }
            return r;
        }
 
        inline const Vec<T,N> unit() const {
            T n = mag();
            Vec<T,N> r;
            for (index_t i = 0; i < N; i++){
                r[i] = this->get(i) / n;
            }
            return r;
        }

        inline T dot(const Vec<T,N> &v) const {
            T sum = 0;
            for (index_t i = 0; i < N; i++){
                sum += this->get(i) * v[i];
            }
            return sum;
        }
        
        //vector magnitude
        inline T mag() const {
            return sqrt(this->mag2());
        }
        
        //magnitude squared
        inline T mag2() const {
            T sum = 0;
            for (index_t i = 0; i < N; i++){
                T cur = this->get(i);
                sum += cur*cur;
            }
            return sum;
        }

        T dist(const Vec<T,N> &pt) const {
            return (this->minus(pt)).mag();
        }

        T dist2(const Vec<T,N> &pt) const {
            return (this->minus(pt)).mag2();
        }

        const Vec<T,N> reflect(const Vec<T,N> &axis) const {
            axis = axis.unit();
            return this->neg() + axis * 2 * this->dot(axis);
        }

        const Vec<T,N> bounce(const Vec<T,N> &normal) const {
            return -reflect(normal);
        }
        
        const Vec<T,N> projectOn(const Vec<T,N> &v) const {
            return v * this->dot(v) / v.mag2();
        }

        const Vec<T,N> mix(const Vec<T,N> &v, T mix) const {
            T xim = 1 - mix;
            Vec<T,N> r;
            for (index_t i = 0; i < N; i++){
                r[i] = xim*this->get(i) + mix*v[i];
            }
            return r;
        }

        T angleTo(const Vec<T,N> &vv) const {
            return acos(vv.unit().dot(unit()));
        }
        
        ///////////////////////
        // Clamping/Rounding //
        ///////////////////////
        
        const Vec<T,N> abs() const {
            Vec<T,N> r;
            for (index_t i = 0; i < N; i++){
                T coord = this->get(i);
                r[i] = coord < 0 ? -coord : coord;
            }
            return r;
        }
        
        //todo: handle integer cases
        const Vec<T,N> floor() const {
            Vec<T,N> r;
            for (index_t i = 0; i < N; i++){
                r[i] = std::floor(this->get(i));
            }
            return r;
        }
        
        //todo: handle integer cases
        const Vec<T,N> ceil() const {
            Vec<T,N> r;
            for (index_t i = 0; i < N; i++){
                r[i] = std::ceil(this->get(i));
            }
            return r;
        }

        const Vec<T,N> min(const Vec<T,N> &v) const {
            Vec<T,N> r;
            for (index_t i = 0; i < N; i++){
                r[i] = std::min(this->get(i), v[i]);
            }
            return r;
        }

        const Vec<T,N> max(const Vec<T,N> &v) const {
            Vec<T,N> r;
            for (index_t i = 0; i < N; i++){
                r[i] = std::max(this->get(i), v[i]);
            }
            return r;
        }
        
        const Vec<T,N> clamp(const Vec<T,N> &one, const Vec<T,N> &two) const {
            Vec<T,N> r;
            for (index_t i = 0; i < N; i++){
                T low = std::min(one[i], two[i]);
                T hi  = std::max(one[i], two[i]);
                r[i]  = std::min(std::max(this->get(i), low), hi);
            }
            return r;
        }
        
        const Vec<T,N> round() const {
            Vec<T,N> r;
            for (index_t i = 0; i < N; i++){
                r[i] = std::floor(this->get(i) + 0.5);
            }
            return r;
        }
        
        template <index_t M> inline Vec<T,M> resized() const {
            Vec<T,M> out;
            const T *p = this->begin();
            std::copy(p, p + std::min(M,N), out.begin());
            return out;
        }

        bool isZero() const {
            for (index_t i = 0; i < N; i++){
                if (this->get(i) != (T)0) return false;
            }
            return true;
        }

        inline index_t getSize() const {
            return N;
        }
        
        //TODO: better hash function. 
        //TODO: doesn't work with large T.
        size_t hashcode() const {
            size_t h = 0;
            const size_t magic = 31;
            const int elembits  = sizeof(T)*8; 
            const int sizetbits = std::numeric_limits<size_t>::digits;
            const int foldings = (elembits + sizetbits - 1) / sizetbits; //ceil(elembits/sizetbits)
            
            typedef typename boost::uint_t<elembits>::least container_t;
            
            for (index_t axis = 0; axis < N; axis++){
                //treat the bits of each component as an unsigned int of sufficient size
                container_t xx = *((container_t*)&(this->get(axis)));
                //take an N-bit float/int and fold it into a possibly 32-bit (or even 16-bit?) size_t
                size_t folded_xx = 0;
                if (elembits > sizetbits){
                    for (index_t i = 0; i < foldings; i++){
                        folded_xx = folded_xx ^ (xx & std::numeric_limits<size_t>::max());
                        xx = xx >> sizetbits;
                    }
                } else {
                    folded_xx = xx;
                }
                h = h*magic + folded_xx;
            }
            return h;
        }

    }; /* end VecCommon class */

    template <typename T, index_t N> const Vec<T,N> VecCommon<T,N>::ones  = Vec<T,N>((T)1);
    template <typename T, index_t N> const Vec<T,N> VecCommon<T,N>::zeros = Vec<T,N>((T)0);

    // used by fromRGB()
    // for integer vectors, you don't want to divide by 255, essentially
    // these perform that logic switch:
    
    template <typename T, typename Enable=void>
    struct _RGBChannelConversion {
        const static double scale = 1 / 255.0;
    };
    
    template <typename T>
    struct _RGBChannelConversion<T, typename boost::enable_if<boost::is_integral<T>,void>::type> {
        const static int scale = 1;
    };
    
}; /* end namespace detail */

} /* end namepsace geom */


/*********************************
 * Hashing                       *
 *********************************/

//allows Vec<T,N>s to work with std::tr1 hash containers
namespace std { namespace tr1 {

   template <typename T, index_t N> struct hash< geom::Vec<T,N> > : public unary_function<geom::Vec<T,N>, size_t> {
       size_t operator()(const geom::Vec<T,N>& v) const {
           return v.hashcode();
       }
   };
   
}}

#endif /* VECBASE_H_ */
