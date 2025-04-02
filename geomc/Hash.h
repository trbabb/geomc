#pragma once

#include <string_view>
#include <functional>
#include <typeindex>

#include <geomc/function/Utils.h>

namespace geom {

template <typename T, typename H>
H hash(const T& obj);


namespace detail {

void siphash(
        const void*    in_bytes,
        const size_t   in_len,
        const void*    key,
              uint8_t* out,
        const size_t   out_len);

} // namespace detail


template <std::integral H>
constexpr H truncated_constant(uint64_t k0, uint64_t k1) {
    static_assert(sizeof(H) <= 16, "truncated_constant() only supports up to 128 bits");
    if constexpr (sizeof(H) > 8) {
        return ((H) k0) << 64 | ((H) k1);
    } else {
        return (H) k1;
    }
}


template <typename T, typename H>
constexpr H type_constant() {
    size_t t_id = std::type_index(typeid(T)).hash_code();
    return truncated_constant<H>(0xdb8f3e91a1d0cde7, t_id);
}

/**
 * @brief Combine two hashes into one.
 */
template <std::unsigned_integral H>
inline constexpr H hash_combine(H h0, H h1) {
    // todo: strongly consider using rabin-karp here.
    //   - multipliers can be taken from "good" lcgs
    //   - for 128 bit, you can use from PCG:
    //     0x2360ed051fc65da4'4385df649fccf645
    if constexpr (sizeof(H) == 16) {
        // constant 0x9e3779b97f4a7c15f39cc0605cedc834 calculated
        // with a multiprecision library; k = 2^128 / phi
        constexpr H k_hi = 0x9e3779b97f4a7c15ULL;
        constexpr H k_lo = 0xf39cc0605cedc834ULL;
        constexpr H k = (k_hi << 64) | k_lo;
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
            Hs >= 2 and Hs <= 16 and is_power_of_two(Hs),
            "bit width not supported for hashing"
        );
    }
}

template <typename H>
inline H hash_bytes(H nonce, const void* obj, size_t size) {
    const char* bytes = reinterpret_cast<const char*>(obj);
    if constexpr (sizeof(size_t) >= sizeof(H)) {
        std::string_view s{bytes, size};
        return std::hash<std::string_view>{}(s) ^ nonce;
    } else {
        uint64_t key[2];
        uint64_t out[2];
        key[0] = nonce;
        if constexpr (sizeof(H) > 8) {
            key[1] = nonce >> 64;
        }
        detail::siphash(
            obj,
            size,
            key,
            reinterpret_cast<uint8_t*>(out),
            sizeof(H) <= 8 ? 8 : 16
        );
        if constexpr (sizeof(H) <= 8) {
            return out[0];
        } else {
            return (((H) out[0]) << 64) | out[1];
        }
    }
}

template <typename H>
inline constexpr H hash_combine_many(H h) {
    return h;
}

template <typename H, typename... Hs>
inline constexpr H hash_combine_many(H h, Hs... hashes) {
    constexpr size_t N = sizeof...(Hs);
    H hs[N] {hashes...};
    return hash_bytes<H>(h, hs, N * sizeof(H));
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


/**
 * @brief Specialize for `T` to disable the use of `std::hash` for `Digest<T>`.
 *
 * This can be used to disambiguate specializations of `Digest` with certain
 * templated custom types. (Most templates will not be ambiguous).
 */
template <typename T>
struct disable_std_hash_fallback : std::false_type {};


// if std::hash is defined for T and there are enough
// bits in size_t to fill a H, use that
template <typename T, typename H>
requires (
    sizeof(H) <= sizeof(size_t)
    and not disable_std_hash_fallback<T>::value
)
struct Digest<T, H> {
    H operator()(const T& obj) const {
        return std::hash<T>{}(obj);
    }
};


// if T literally H, use identity
template <typename H>
struct Digest<H,H> {
    H operator()(const H& h) const {
        return h;
    }
};


// digest an arithmetic type by casting its bits.
// if H is larger than size_t, we can't use std::hash,
// but if T is bigger than size_t and smaller than H,
// we'll just use its bits directly.
template <typename T, typename H>
requires (
        sizeof(H) > sizeof(size_t)  // std::hash cannot fill an H
    and sizeof(H) > sizeof(T)       // not enough bits in T to fill an H,
    and sizeof(T) > sizeof(size_t)  // but T is bigger than size_t,
                                    // so we get more bits than std::hash;
                                    // might as well use all available bits directly
    and (std::is_arithmetic_v<H> or std::is_enum_v<H>)
)
struct Digest<T, H> {
    H operator()(const T& x) const {
        // todo: do we want to mix up the bits rather than directly using them?
        //   note that we might want, e.g., all ints to hash to the same thing,
        //   because promotion could lead to confusing results otherwise.
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


// digest an arithmetic type by folding its bits.
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
    H operator()(T x) const {
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


// digest a string
template <typename H>
requires (sizeof(H) > sizeof(size_t))
struct Digest<std::string, H> {
    H operator()(const std::string& s) const {
        // should match std::string view below:
        H nonce = truncated_constant<H>(0xf10749c80725844e, 0xe88cfd8ec23556e4);
        return hash_bytes<H>(0, s.data(), s.size());
    }
};


// digest a string view
template <typename H>
requires (sizeof(H) > sizeof(size_t))
struct Digest<std::string_view, H> {
    H operator()(const std::string_view& s) const {
        // should match std::string above:
        H nonce = truncated_constant<H>(0xf10749c80725844e, 0xe88cfd8ec23556e4);
        return hash_bytes<H>(0, s.data(), s.size());
    }
};


// digest a tuple
template <typename... Ts, typename H>
struct Digest<std::tuple<Ts...>, H> {
    H operator()(const std::tuple<Ts...>& t) const {
        return std::apply(
            [](const Ts&... elems) {
                return hash_combine_many(
                    truncated_constant<H>(0xf804bbd33c7e7769, 0x439588acee006cd0),
                    geom::hash<Ts, H>(elems)...
                );
            },
            t
        );
    }
};


// digest an optional
template <typename T, typename H>
struct Digest<std::optional<T>, H> {
    H operator()(const std::optional<T>& opt) const {
        // hash the contained value with a nonce unique to optional;
        // or a nonce unique to the contained type if empty
        return geom::hash_combine<H>(
            geom::truncated_constant<H>(0x40949d0dce18bcac, 0xc68f7668569db4f4),
            opt ? geom::hash<T,H>(*opt) : geom::type_constant<T,H>()
        );
    }
};


// digenst a variant
template <typename... Ts, typename H>
struct Digest<std::variant<Ts...>, H> {
    H operator()(const std::variant<Ts...>& v) const {
        return geom::hash_combine_many<H>(
            // nonce unique to variant
            geom::truncated_constant<H>(0xc073bebcbe0daa53, 0x9100685c51aa2de4),
            // entropy unique to the active alternative
            geom::hash<H>(v.index()),
            // hash the value
            std::visit(
                [](const auto& x) {
                    return geom::hash(x);
                },
                v
            )
        );
    }
};


/// Produce a hash of an object.
template <typename T, typename H>
H hash(const T& obj) {
    return Digest<T,H>{}(obj);
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


} // namespace geom
