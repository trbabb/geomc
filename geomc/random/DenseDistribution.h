#pragma once

#include <random>
#include <geomc/shape/Rect.h>

namespace geom {

/**
 * @addtogroup random
 * @{
 */

/**
 * @brief A random number generator that produces uniformly-distributed values.
 *
 * Many uniform random number generators only generate a small subset of the
 * representable floating point values. This generator produces all possible
 * values in the range [0,1] with equal probability by directly constructing
 * the mantissa and exponent of a floating point number.
 *
 * Pharr and Humphreys discuss the method used here:
 *   https://pharr.org/matt/blog/2022/03/05/sampling-fp-unit-interval
 *
 * Non-float types fall back to std::uniform_real_distribution.
 */
template <typename T>
struct DenseUniformDistribution : public std::uniform_real_distribution<T> {};

/// @brief Dense uniform distribution specialization for float.
template <>
struct DenseUniformDistribution<float> {
private:
    Rect<float,1> _range;
    
public:
    using param_type  = Rect<float,1>;
    using result_type = float;
    
    constexpr DenseUniformDistribution():
        _range(0.f, 1.f) {}
    constexpr DenseUniformDistribution(Rect<float,1> range):
        _range(range) {}
    constexpr DenseUniformDistribution(float a, float b):
        _range(Rect<float,1>::spanning_corners(a, b)) {}
    
    template <typename Generator>
    float operator()(Generator& g) const {
        // without IEEE754, our bit twiddling will not work.
        // generally this will only fail for "exotic" systems:
        static_assert(std::numeric_limits<float>::is_iec559, "floating point must be ieee754");
        
        union FPBox {
            float f;
            uint32_t i;
        };
        
        uint32_t mant;
        int exp, hi_exp, lo_exp;
        FPBox lo, hi, ans;
        
        lo.f = 0.0;
        hi.f = 1.0;
        
        lo_exp = (lo.i >> 23) & 0xFF;
        hi_exp = (hi.i >> 23) & 0xFF;
        
        // xxx: use clz()!
        auto randbit = std::uniform_int_distribution<uint32_t>();
        uint32_t r = randbit(g);
        //not >= because exp is decremented at the end of the last loop
        for (exp = hi_exp - 1; exp > lo_exp; exp--) {
            if (r & 1) break;
            r >>= 1;
            // we could decrement up to 127 times, but we only have 32 bits of
            // randomness at a time. generate more bits.
            if (exp == hi_exp - 32) [[unlikely]] {
                r = randbit(g);
            }
        }
        
        mant = (randbit(g) & 0xFFFFFE00) >> 9; // use the high quality bits
        
        if (mant == 0 and (r & 1)) exp++; // border values must not be skewed towards 1
        
        ans.i = (((uint32_t)exp) << 23) | mant; // combine exp and mantissa
        return _range.remap(ans.f);
    }
    
    void reset() {}
    Rect<float,1> param() const { return _range; }
    void param(Rect<float,1> range) { _range = range; }
    
    float min() const { return _range.lo; }
    float max() const { return _range.hi; }
    
    float a() const { return _range.lo; }
    float b() const { return _range.hi; }
    
    bool operator==(const DenseUniformDistribution& other) const {
        return _range == other._range;
    }
};

/// @brief Dense uniform distribution specialization for double.
template <>
struct DenseUniformDistribution<double> {
private:
    Rect<double,1> _range;
    
public:
    using param_type  = Rect<double,1>;
    using result_type = double;
    
    constexpr DenseUniformDistribution():
        _range(0., 1.) {}
    constexpr DenseUniformDistribution(Rect<double,1> range):
        _range(range) {}
    constexpr DenseUniformDistribution(double a, double b):
        _range(Rect<double,1>::spanning_corners(a, b)) {}
    
    template <typename Generator>
    double operator()(Generator& g) const {
        union FPBox {
            double   d;
            uint64_t l;
        };
        
        // without IEEE754, our bit twiddling will not work.
        // generally this will only fail for "exotic" systems:
        static_assert(std::numeric_limits<double>::is_iec559, "floating point must be ieee754");
        
        int exp, hi_exp, lo_exp;
        uint64_t mant;
        FPBox lo, hi, ans;
        
        lo.d = 0.0;
        hi.d = 1.0;
        
        lo_exp = (lo.l >> 52) & 0x7FF;
        hi_exp = (hi.l >> 52) & 0x7FF;
        
        auto randbit = std::uniform_int_distribution<uint64_t>();
        uint64_t r = randbit(g);
        //not >= because exp is decremented at the end of the last loop
        for (exp = hi_exp - 1; exp > lo_exp; exp--) {
            if (r & 1) break;
            r >>= 1;
            // we could decrement up to 1023 times, but we only have 64 bits of
            // randomness at a time. generate more bits. this is EXTREMELY unlikely.
            if (exp == hi_exp - 64) [[unlikely]] {
                r = randbit(g);
            }
        }
        
        mant = (randbit(g) & 0xFFFFFFFFFFFFF000LL) >> 12; // use the high quality bits
        
        if (mant == 0 && (r & 1)) exp++; // border values must not be skewed towards 1
        
        ans.l = (((uint64_t)exp) << 52) | mant; // combine exp and mantissa
        return _range.remap(ans.d);
    }
    
    void reset() {}
    Rect<double,1> param() const { return _range; }
    void param(Rect<double,1> range) { _range = range; }
    
    double min() const { return _range.lo; }
    double max() const { return _range.hi; }
    
    double a() const { return _range.lo; }
    double b() const { return _range.hi; }
    
    bool operator==(const DenseUniformDistribution& other) const {
        return _range == other._range;
    }
};


/// @brief Dense uniform distribution specialization for Duals.
template <typename T, DiscontinuityPolicy P>
struct DenseUniformDistribution<Dual<T,P>> {
private:
    DenseUniformDistribution<T> _d;
    Rect<Dual<T,P>,1> _range;
public:
    using dual_t      = Dual<T,P>;
    using param_type  = Rect<dual_t,1>;
    using result_type = dual_t;
    
    constexpr DenseUniformDistribution():
        _range(0., 1.) {}
    constexpr DenseUniformDistribution(param_type range):
        _range(range) {}
    constexpr DenseUniformDistribution(T a, T b):
        _range(Rect<dual_t,1>::spanning_corners(a, b)) {}
    
    template <typename Generator>
    Dual<T,P> operator()(Generator& g) const {
        dual_t d = _d(g);
        return _range.remap(d);
    }
    
    void reset() { _d.reset(); }
    param_type param() const { return _range; }
    void param(param_type range) { _range = range; }
    
    dual_t min() const { return _range.lo; }
    dual_t max() const { return _range.hi; }
    
    T a() const { return _range.lo.real(); }
    T b() const { return _range.hi.real(); }
    
    bool operator==(const DenseUniformDistribution& other) const {
        return _range == other._range;
    }
};

/// @}

} // namespace geom
