#pragma once

#include <random>
#include <bit>
#include <geomc/shape/Rect.h>

// todo: we should use high order bits (higher quality) for the exponent
//    and low order bits for the mantissa. exponent is more important for
//    the distribution of the result.

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
 * Pharr discusses the method used here:
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
        _range(Rect<float,1>::from_corners(a, b)) {}
    
    template <typename Generator>
    float operator()(Generator& g) const {
        union FPBox {
            float    f;
            uint32_t i;
        };
        
        // without IEEE754, our bit twiddling will not work.
        // generally this will only fail for "exotic" systems:
        static_assert(std::numeric_limits<float>::is_iec559, "floating point must be ieee754");
        
        auto randbits = std::uniform_int_distribution<uint64_t>();
        constexpr size_t mantissa_bits = 23;
        size_t scale_bits = 64 - mantissa_bits;
        uint64_t r = randbits(g);
        int32_t  exp = 126;
        uint32_t mant = r >> scale_bits;
        uint32_t z;
        while (true) {
            z   = std::countr_zero(r);
            exp = std::max<int32_t>(exp - z, 0);
            if (z < scale_bits or exp == 0) [[likely]] break;
            // we will get here only with probability 2^-40
            scale_bits = 64;
            r = randbits(g);
        }
        if (mant == 0) [[unlikely]] {
            // avoid biasing the result toward zero. above, it is not
            // possible to generate 1.0. with probability 1/2, bump the exponent.
            exp += ((r >> (scale_bits - 1)) & 1);
        }
        FPBox ans;
        ans.i = (exp << mantissa_bits) | mant; // combine exp and mantissa
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
        _range(Rect<double,1>::from_corners(a, b)) {}
    
    template <typename Generator>
    double operator()(Generator& g) const {
        union FPBox {
            double   d;
            uint64_t l;
        };
        
        // without IEEE754, our bit twiddling will not work.
        // generally this will only fail for "exotic" systems:
        static_assert(std::numeric_limits<float>::is_iec559, "floating point must be ieee754");
        
        auto randbits = std::uniform_int_distribution<uint64_t>();
        constexpr size_t mantissa_bits = 52;
        size_t scale_bits = 64 - mantissa_bits;
        uint64_t r = randbits(g);
        int32_t  exp = 1022;
        uint64_t mant = r >> scale_bits;
        uint32_t z;
        while (true) {
            z   = std::countr_zero(r);
            exp = std::max<uint32_t>(exp - z, 0);
            if (z < scale_bits or exp == 0) [[likely]] break;
            // we get here with probablity 2^-11
            scale_bits = 64;
            r = randbits(g);
        }
        if (mant == 0) [[unlikely]] {
            // avoid biasing the result toward zero. above, it is not
            // possible to generate 1.0. with probability 1/2, bump the exponent.
            exp += ((r >> (scale_bits - 1)) & 1);
        }
        FPBox ans;
        ans.l = (((uint64_t)exp) << mantissa_bits) | mant; // combine exp and mantissa
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
        _range(Rect<dual_t,1>::from_corners(a, b)) {}
    
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
