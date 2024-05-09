#pragma once

#include <string_view>
#include <functional>

#include <geomc/function/Utils.h>

namespace geom {


template <std::integral H>
constexpr H truncated_constant(uint64_t k0, uint64_t k1) {
    static_assert(sizeof(H) <= 16, "truncated_constant() only supports up to 128 bits");
    if constexpr (sizeof(H) > 8) {
        return ((H) k0) << 64 | ((H) k1);
    } else {
        return (H) k0;
    }
}


/**
 * @brief Partially-specializable hash function object for arbitrary
 * object and hash types.
 * 
 * To enable hashing of a custom type, specialize this template for that type.
 * In general it is good to hash with a type-specific nonce so that distinct types
 * that have the same member data do not collide.
 */
template <typename T, typename H=std::size_t>
struct Digest {};


// if std::hash is defined for T and there are enough
// bits in size_t to fill a H, use that
template <typename T, typename H>
requires (sizeof(H) <= sizeof(size_t))
struct Digest<T, H> {
    H operator()(const T& obj) const {
        return std::hash<T>{}(obj);
    }
};


// if H is larger than size_t, we can't use std::hash,
// but if T is bigger than size_t and smaller than H,
// we'll just use its bits directly.
template <typename T, typename H>
requires (
        sizeof(H) > sizeof(size_t)  // std::hash cannot fill an H
    and sizeof(H) > sizeof(T)       // not enough bits in T to fill an H
    and sizeof(T) > sizeof(size_t)  // but T is bigger than size_t,
                                    // so we get more bits than std::hash;
                                    // might as well use all available bits directly
    and (std::is_arithmetic_v<H> or std::is_enum_v<H>)
)
struct Digest<T, H> {
    H operator()(const T& x) const {
        // todo: do we want to mix up the bits rather than directly using them?
        if constexpr (std::floating_point<T>) {
            x = x + (T)0; // force -0 to +0
            if constexpr (sizeof(T) == 2) {
                return *reinterpret_cast<const uint16_t*>(&x);
            } else if constexpr (sizeof(T) == 4) {
                return *reinterpret_cast<const uint32_t*>(&x);
            } else if constexpr (sizeof(T) == 8) {
                return *reinterpret_cast<const uint64_t*>(&x);
            } else {
                H h = 0;
                const uint8_t* src_bytes = reinterpret_cast<const uint8_t*>(&x);
                      uint8_t* dst_bytes = reinterpret_cast<      uint8_t*>(&h);
                std::copy(src_bytes, src_bytes + sizeof(T), dst_bytes);
                return h;
            }
        } else {
            return (H) x;
        }
    }
};


// if H is larger than size_t, we can't use std::hash; and
// if T is arithmetic and (over-) fills an H, then we hash it by treating it
// as an array of Hs.
template <typename T, typename H>
requires (
    sizeof(H) > sizeof(size_t)  // std::hash cannot fill an H
    and sizeof(H) <= sizeof(T)  // T's bits can fill or overfill an H
    and std::is_arithmetic_v<H> // T not a class type
)
struct Digest<T, H> {
    H operator()(const T& x) const {
        constexpr size_t Ts = sizeof(T);
        constexpr size_t Hs = sizeof(H);
        // this ought to always be possible because arithmetic types
        // are generally always a power of two number of bytes
        static_assert(
            Ts % Hs == 0,
            "hash_combine expects that sizeof(T) is a multiple of sizeof(H)"
        );
        if constexpr (std::floating_point<T>) {
            x = x + (T)0; // force -0 to +0
        }
        const H* v = reinterpret_cast<const H*>(&x);
        H h = v[0];
        for (size_t i = 1; i < Ts / Hs; ++i) {
            h = hash_combine(h, v[i]);
        }
        return h;
    }
};


/// Produce a hash of an object.
template <typename T, typename H>
H hash(const T& obj) {
    return Digest<T,H>{}(obj);
}

/**
 * @brief Combine two hashes into one.
 */
template <std::unsigned_integral H>
inline constexpr H hash_combine(H h0, H h1) {
    if constexpr (sizeof(H) == 16) {
        // constant 0x9e3779b97f4a7c15f39cc0605cedc834 calculated
        // with a multiprecision library; k = 2^128 / phi
        constexpr H k_hi = 0x9e3779b97f4a7c15ULL;
        constexpr H k_lo = 0xf39cc0605cedc834ULL;
        H k = (k_hi << 64) | k_lo;
        h0 ^= h1 + k + (h0 << 24) + (h0 >> 8);
        return h0;
    } else if constexpr (sizeof(H) == 8) {
        h0 ^= h1 + 0x9e3779b97f4a7c15ULL + (h0 << 12) + (h0 >> 4);
        return h0;
    } else if constexpr (sizeof(H) == 4) {
        // boost's hash_combine. everyone uses this, though there are reports that it's shit.
        // see https://stackoverflow.com/questions/8513911/,
        //     https://stackoverflow.com/questions/5889238/
        // i am not a cryptographer, though, so I don't expect to do much better myself.
        // we should make this better, though
        h0 ^= h1 + 0x9e3779b9U + (h0 << 6) + (h0 >> 2);
        return h0;
    } else if constexpr (sizeof(H) == 2) {
        // same as above. ugh
        h0 ^= h1 + 0x9e37U + (h0 << 3) + (h0 >> 1);
        return h0;
    } else {
        constexpr size_t Hs = sizeof(H);
        static_assert(
            Hs >= 2 and Hs <= 8 and is_power_of_two(Hs),
            "bit width not supported for hashing"
        );
    }
}

template <typename H>
inline constexpr H hash_combine_many(H h) {
    return h;
}

template <typename H, typename... Hs>
inline constexpr H hash_combine_many(H h, Hs... hashes) {
    return hash_combine(h, hash_combine_many(hashes...));
}

template <typename T, typename H>
inline H hash_array(H nonce, const T* objs, size_t count) {
    H h = nonce;
    for (size_t i = 0; i < count; ++i) {
        h = hash_combine(h, geom::hash<T,H>(objs[i]));
    }
    return h;
}

template <typename H, typename... Ts>
inline H hash_many(H nonce, const Ts&... objs) {
    return hash_combine_many(nonce, geom::hash<Ts,H>(objs)...);
}

/*
template <typename H>
inline size_t hash_bytes(const void* obj, size_t size) {
    const char* bytes = reinterpret_cast<const char*>(obj);
    if constexpr (sizeof(size_t) >= sizeof(H)) {
        std::string_view s{bytes, size};
        return std::hash<std::string_view>{}(s);
    } else {
        // xxx: need to hash arbitrary bytes to an H
    }
}
*/


} // namespace geom
