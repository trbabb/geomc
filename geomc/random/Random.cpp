/*
 * Random.cpp
 *
 *  Created on: Feb 22, 2009
 *      Author: Tim Babb
 */

#include <boost/static_assert.hpp>
#include <assert.h>
#include <algorithm>
#include <limits>

#include <geomc/random/Random.h>

#define INT_SIGNBIT       (1 << std::numeric_limits<int>::digits)
#define LONG_SIGNBIT      (1L << std::numeric_limits<long>::digits)
#define LONG_LONG_SIGNBIT (1LL << std::numeric_limits<long long>::digits)

using namespace std;
using namespace geom;

typedef union fpbox {
    uint32_t i;
    float f;
    uint64_t l;
    double d;
} FPBox;

/**********************************
 * Structors                      *
 **********************************/

Random::Random():_bitpool(0),
                 _bitsleft(0){}

Random::~Random(){}

/**********************************
 * Elementary Specializations     *
 **********************************/

namespace geom {

///@ingroup random 
///@{

/// @return A random boolean value.
template <> bool Random::rand<bool>(){
    if (_bitsleft == 0){
        _bitpool = this->rand32();
        _bitsleft = 32;
    }
    int bit = _bitpool & 1;
    _bitpool = _bitpool >> 1;
    _bitsleft--;
    return bit;
}

/**
 * @return An `unsigned int` between 0 and `UINT_MAX`.
 */
template <> unsigned int Random::rand<unsigned int>(){
    unsigned int bits = 0;
    const int chunks = (std::numeric_limits<unsigned int>::digits + 32 - 1) / 32; //ceiling(x/y) = (x + y - 1) / y
#if UINT_MAX > (1LL << 32) - 1
    for (int i = 0; i < chunks; i++){
        bits = (bits << 32) | rand32();
    }
#else
    bits =(unsigned int)rand32();
#endif
    return bits;
}

/**
 * @return An `int` between `INT_MIN` and `INT_MAX`.
 */
template <> int Random::rand<int>(){
    return rand<unsigned int>(); //MSB is now a sign bit; covers the entire range of <int> type. might we want to change this? 
}

/**
 * @return An `unsigned long` between 0 and `ULONG_MAX`.
 */
template <> unsigned long Random::rand<unsigned long>(){
    unsigned long bits = 0;
    const int chunks = (std::numeric_limits<unsigned long>::digits + 32 - 1) / 32; //ceiling(x/y) = (x + y - 1) / y
#if ULONG_MAX > (1LL << 32) - 1
    for (int i = 0; i < chunks; i++){
        bits = (bits << 32) | rand32();
    }
#else
    bits = (unsigned long)rand32();
#endif
    return bits;
}

/**
 * @return A `long` between `LONG_MIN` and `LONG_MAX`.
 */
template <> long Random::rand<long>(){
    return rand<unsigned long>();
}

/**
 * @return An `unsigned long long` betwen 0 and `ULLONG_MAX`.
 */
template <> unsigned long long Random::rand<unsigned long long>(){
    unsigned long long bits = 0;
    const int chunks = (std::numeric_limits<unsigned long long>::digits + 32 - 1) / 32; //ceiling(x/y) = (x + y - 1) / y
    for (int i = 0; i < chunks; i++){
        bits = (bits << 32) | rand32();
    }
    return bits;
}

/**
 * @return A `long long` between `LLONG_MIN` and `LLONG_MAX`.
 */
template <> long long Random::rand<long long>(){
    return rand<unsigned long long>();
}

/* From the unpublished paper by Allen B. Downey.
 *
 * The conventional method of dividing a random number
 * by a constant excludes about 93% of the representable
 * floating point values between 0.0 and 1.0. Instead
 * we pick the bits of our number explicitly, but must
 * be careful to do so such that the distribution is
 * uniform. This method is fast, and makes efficient
 * use of our pseudorandom bits.
 *
 * The gist: Pick the mantissa, then pick the exponent
 * by iterative coin-flipping.
 */

//TODO: use numeric_limits instead of hard-coded masks
//see: boost::integer's int_t<bits>

/**
 * Generate a floating-point value with uniform distribution between 0 and 1.0.
 * 
 * This method is preferable to use over the common practice of 
 * `rand<int>() / (float)INT_MAX`, as the latter excludes a large fraction of 
 * the representable floating point numbers between 0.0 and 1.0 (particularly those
 * near zero), resulting in reduced entropy. 
 * 
 * This method produces a number by carefully picking the bits of the exponent 
 * and mantissa explicitly, in such a way that the distibution remains uniform.
 */
template <> float Random::rand<float>(){
    // without IEEE754, our bit twiddling will not work.
    // generally this will only fail for "exotic" systems:
    BOOST_STATIC_ASSERT(std::numeric_limits<float>::is_iec559);
    
    uint32_t mant;
    int exp, hi_exp, lo_exp;
    FPBox lo, hi, ans;

    lo.f = 0.0;
    hi.f = 1.0;

    lo_exp = (lo.i >> 23) & 0xFF;
    hi_exp = (hi.i >> 23) & 0xFF;

    //not >= because exp is decremented at the end of the last loop
    for (exp = hi_exp - 1; exp > lo_exp; exp--){
        if (rand<bool>()) break;
    }

    mant = (rand32() & 0xFFFFFE00) >> 9; //use the high quality bits

    if (mant == 0 && rand<bool>()) exp++; //border values must not be skewed towards 1

    ans.i = (((uint32_t)exp) << 23) | mant; //combine exp and mantissa
    return ans.f;
}

/**
 * Generate a floating-point value with uniform distribution between 0 and 1.0.
 * 
 * This method is preferable to use over the common practice of 
 * `rand<int>() / (float)INT_MAX`, as the latter excludes a large fraction of 
 * the representable floating point numbers between 0.0 and 1.0 (particularly those
 * near zero), resulting in reduced entropy. 
 * 
 * This method produces a number by carefully picking the bits of the exponent 
 * and mantissa explicitly, in such a way that the distibution remains uniform.
 */
template <> double Random::rand<double>(){
    // without IEEE754, our bit twiddling will not work.
    // generally this will only fail for "exotic" systems:
    BOOST_STATIC_ASSERT(std::numeric_limits<double>::is_iec559);
    
    int exp, hi_exp, lo_exp;
    uint64_t mant;
    FPBox lo, hi, ans;

    lo.d = 0.0;
    hi.d = 1.0;

    lo_exp = (lo.l >> 52) & 0x7FF;
    hi_exp = (hi.l >> 52) & 0x7FF;

    //not >= because exp is decremented at the end of the last loop
    for (exp = hi_exp - 1; exp > lo_exp; exp--){
        if (rand<bool>()) break;
    }
    
    mant = ((uint64_t)rand32() << 32) | rand32();
    mant = (mant & 0xFFFFFFFFFFFFF000LL) >> 12; //use the high quality bits

    if (mant == 0 && rand<bool>()) exp++; //border values must not be skewed towards 1

    ans.l = (((uint64_t)exp) << 52) | mant; //combine exp and mantissa
    return ans.d;
}

/**********************************
 * Positive Range Specializations *
 **********************************/

/**
 * @param hi Upper bound of random number.
 * @return An `unsigned int` between 0 and `hi`.
 */
template <> unsigned int Random::rand<unsigned int>(unsigned int hi){
    unsigned int bits, val, lo;

    do {
        bits = rand<unsigned int>();
        val = bits % hi;
        lo = bits - val;
        //while this chunk wraps around to zero:
        //i.e. lo sign bit is 1, end-of-chunk sign bit is 0
    } while ( (lo & INT_SIGNBIT) & ~(lo+(hi-1)) );
    return val;
}

/**
 * @param hi Upper or lower bound of random number.
 * @return A number between 0 and `hi` (regardless of the sign of `hi`).
 */
template <> int Random::rand<int>(int hi){
    if (hi < 0){
        return -rand<unsigned int>(-hi);
    } else {
        return rand<unsigned int>(hi);
    }
}

/**
 * @param hi Upper bound of random number.
 * @return An `unsigned long` between 0 and `hi`.
 */
template <> unsigned long Random::rand<unsigned long>(unsigned long hi){
    unsigned long bits, val, lo;

    do {
        bits = rand<unsigned long>();
        val = bits % hi;
        lo = bits - val;
        //while this chunk wraps around to zero:
        //i.e. lo sign bit is 1, end-of-chunk sign bit is 0
    } while ( (lo & LONG_SIGNBIT) & ~(lo+(hi-1)) );
    return val;
}

/**
 * @param hi Upper or lower bound of random number.
 * @return A `long` between 0 and `hi`, regardless of the sign of `hi`.
 */
template <> long Random::rand<long>(long hi){
    if (hi < 0){
        return -rand<unsigned long>(-hi);
    } else {
        return rand<unsigned long>(hi);
    }
}

/**
 * @param hi Upper bound of random number.
 * @return An `unsigned long long` between 0 and `hi`.
 */
template <> unsigned long long Random::rand<unsigned long long>(unsigned long long hi){
    unsigned long long bits, val, lo;

    do {
        bits = rand<unsigned long long>();
        val = bits % hi;
        lo = bits - val;
        //while this chunk wraps around to zero:
        //i.e. lo sign bit is 1, end-of-chunk sign bit is 0
    } while ( (lo & LONG_LONG_SIGNBIT) & ~(lo+(hi-1)) );
    return val;
}

/**
 * @param hi Upper or lower bound of random number.
 * @return A `long long` between 0 and `hi`, regardless of the sign of `hi`.
 */
template <> long long Random::rand<long long>(long long hi){
    if (hi < 0){
        return -rand<unsigned long long>(-hi);
    } else {
        return rand<unsigned long long>(hi);
    }
}

/**
 * @param hi Upper or lower bound of random number.
 * @return A `float` between 0 and `hi`, regardless of the sign of `hi`.
 */
template <> float Random::rand<float>(float hi){
    return hi * rand<float>();
}

/**
 * @param hi Upper or lower bound of random number.
 * @return A `double` between 0 and `hi`, regardless of the sign of `hi`.
 */
template <> double Random::rand<double>(double hi){
    return hi * rand<double>();
}

/**********************************
 * Full Range Specializations     *
 **********************************/

// todo: some of these will fail for ranges crossing zero where lo - hi is
// greater than INT_MAX.

/**
 * @param lo Lower bound of random number.
 * @param hi Upper bound of random number.
 * @return An `int` between `lo` and `hi`.
 */
template <> int Random::rand<int>(int lo, int hi){
    int lo1 = min(hi,lo);
    return rand<int>(max(hi,lo) - lo1) + lo1;
}

/**
 * @param lo Lower bound of random number.
 * @param hi Upper bound of random number.
 * @return A `long long` between `lo` and `hi`.
 */
template <> long long Random::rand<long long>(long long lo, long long hi){
    long long lo1 = min(hi,lo);
    return rand<long long>(max(hi,lo) - lo1) + lo1;
}

/**
 * @param lo Lower bound of random number.
 * @param hi Upper bound of random number.
 * @return An `float` between `lo` and `hi`.
 */
template <> float Random::rand<float>(float lo, float hi){
    return (hi-lo)*rand<float>() + lo;
}

/**
 * @param lo Lower bound of random number.
 * @param hi Upper bound of random number.
 * @return An `double` between `lo` and `hi`.
 */
template <> double Random::rand<double>(double lo, double hi){
    return (hi-lo)*rand<double>() + lo;
}

/// @} //ingroup random

} //end namespace geom
