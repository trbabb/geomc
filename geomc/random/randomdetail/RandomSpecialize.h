/* 
 * File:   RandomSpecialize.h
 * Author: tbabb
 *
 * Created on June 3, 2014, 11:04 PM
 */

#ifndef RANDOMSPECIALIZE_H
#define	RANDOMSPECIALIZE_H

namespace geom {

/// @return A random boolean value.
template <> bool Random::rand<bool>();

/**
 * @return An `unsigned int` between 0 and `UINT_MAX`.
 */
template <> unsigned int Random::rand<unsigned int>();

/**
 * @return An `int` between `INT_MIN` and `INT_MAX`.
 */
template <> int Random::rand<int>();

/**
 * @return An `unsigned long` between 0 and `ULONG_MAX`.
 */
template <> unsigned long Random::rand<unsigned long>();

/**
 * @return A `long` between `LONG_MIN` and `LONG_MAX`.
 */
template <> long Random::rand<long>();

/**
 * @return An `unsigned long long` betwen 0 and `ULLONG_MAX`.
 */
template <> unsigned long long Random::rand<unsigned long long>();

/**
 * @return A `long long` between `LLONG_MIN` and `LLONG_MAX`.
 */
template <> long long Random::rand<long long>();

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
template <> float Random::rand<float>();

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
template <> double Random::rand<double>();

/**
 * @param hi Upper bound of random number.
 * @return An `unsigned int` between 0 and `hi`.
 */
template <> unsigned int Random::rand<unsigned int>(unsigned int hi);

/**
 * @param hi Upper or lower bound of random number.
 * @return A number between 0 and `hi` (regardless of the sign of `hi`).
 */
template <> int Random::rand<int>(int hi);

/**
 * @param hi Upper bound of random number.
 * @return An `unsigned long` between 0 and `hi`.
 */
template <> unsigned long Random::rand<unsigned long>(unsigned long hi);

/**
 * @param hi Upper or lower bound of random number.
 * @return A `long` between 0 and `hi`, regardless of the sign of `hi`.
 */
template <> long Random::rand<long>(long hi);

/**
 * @param hi Upper bound of random number.
 * @return An `unsigned long long` between 0 and `hi`.
 */
template <> unsigned long long Random::rand<unsigned long long>(unsigned long long hi);

/**
 * @param hi Upper or lower bound of random number.
 * @return A `long long` between 0 and `hi`, regardless of the sign of `hi`.
 */
template <> long long Random::rand<long long>(long long hi);

/**
 * @param hi Upper or lower bound of random number.
 * @return A `float` between 0 and `hi`, regardless of the sign of `hi`.
 */
template <> float Random::rand<float>(float hi);

/**
 * @param hi Upper or lower bound of random number.
 * @return A `double` between 0 and `hi`, regardless of the sign of `hi`.
 */
template <> double Random::rand<double>(double hi);

/**
 * @param lo Lower bound of random number.
 * @param hi Upper bound of random number.
 * @return An `int` between `lo` and `hi`.
 */
template <> int Random::rand<int>(int lo, int hi);

/**
 * @param lo Lower bound of random number.
 * @param hi Upper bound of random number.
 * @return A `long long` between `lo` and `hi`.
 */
template <> long long Random::rand<long long>(long long lo, long long hi);

/**
 * @param lo Lower bound of random number.
 * @param hi Upper bound of random number.
 * @return An `float` between `lo` and `hi`.
 */
template <> float Random::rand<float>(float lo, float hi);

/**
 * @param lo Lower bound of random number.
 * @param hi Upper bound of random number.
 * @return An `double` between `lo` and `hi`.
 */
template <> double Random::rand<double>(double lo, double hi);


} // namespace geom

#endif	/* RANDOMSPECIALIZE_H */

