/* 
 * File:   Hash.h
 * Author: tbabb
 * 
 * Code for hashing a string of bytes of arbitrary length to an
 * int of type <size_t>.
 * 
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * Therefore, hashcodes will be system-dependent. In addition to pointer
 * size, endian-ness matters. Do not rely on hashes to be identical
 * across systems!
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *
 * Created on August 3, 2013, 2:06 PM
 */

#ifndef HASH_H
#define	HASH_H

#include <geomc/geomc_defs.h>

namespace geom {
namespace detail {
    
    template <size_t bytes>
    struct HashParams {};

    template <>
    struct HashParams<8> {
        static const size_t multiplier = 6364136223846793005LL;
        static const size_t increment  = 1442695040888963407LL;
    };

    template <>
    struct HashParams<6> {
        static const size_t multiplier = 25214903917;
        static const size_t increment  = 11;
    };

    template <>
    struct HashParams<4> {
        static const size_t multiplier = 1664525;
        static const size_t increment  = 1013904223;
    };

    template <>
    struct HashParams<3> {
        static const size_t multiplier = 1140671485;
        static const size_t increment  = 12820163;
    };

    template <>
    struct HashParams<2> {
        static const size_t multiplier = 15731; // <- a prime. pulled outta mah butt.
        static const size_t increment  = 52475; // <- chosen by a fair dice roll. guaranteed to be random.
    };

    template <>
    struct HashParams<1> {
        static const size_t multiplier = 97;  // another buttsourced prime
        static const size_t increment  = 136; // this seems like a nice number.
    };
}; // namespace detail


size_t general_hash(const void *data, size_t bytes);

}; // namespace geom

#endif	/* HASH_H */

