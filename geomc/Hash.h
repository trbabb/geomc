#pragma once

#include <string_view>
#include <functional>

#include <geomc/geomc_defs.h>

namespace geom {

// convenience function which deduces the type of the hash function
template <typename T>
inline size_t hash(const T& obj) {
    return std::hash<T>{}(obj);
}

inline size_t hash_combine(size_t h0, size_t h1) {
    if constexpr (sizeof(size_t) == 8) {
        /*
        // this offers a solution which doesn't use
        // a variant of the shitty boost combine: 
        // https://stackoverflow.com/questions/8513911/
        // but doesn't justify it.
        constexpr size_t k = 0x9E3779B97F4A7C15ULL;
        size_t a = (h0 ^ h1) * k;
        a ^= (a >> 47);
        size_t b = (h0 ^ a) * k;
        b ^= (b >> 47);
        return b * k;
        */
       
        // constant from: https://stackoverflow.com/questions/5889238/
        h0 ^= h1 + 0x517cc1b727220a95ULL + (h0 << 12) + (h0 >> 4);
        return h0;
    } else if constexpr (sizeof(size_t) == 4) {
        // boost's hash_combine. everyone uses this even though it's known to be shit.
        // i am not a cryptographer, though, so I don't expect to do much better myself.
        // we should make this better, though
        h0 ^= h1 + 0x9e3779b9U + (h0 << 6) + (h0 >> 2);
        return h0;
    } else if constexpr (sizeof(size_t) == 2) {
        // same as above. ugh
        h0 ^= h1 + 0x9e37U + (h0 << 3) + (h0 >> 1);
        return h0;
    } else {
        static_assert(sizeof(size_t) > 2, "size_t size not supported for hashing");
    }
}

template <typename... Hs>
inline size_t hash_combine_many(size_t h, Hs... hashes) {
    return hash_combine(h, hash_combine_many(hashes...));
}

template <typename T>
inline size_t hash_many(const T* objs, size_t count) {
    size_t h = (size_t) 0xa5668ab063399788ULL; // arbitrary nonce
    for (size_t i = 0; i < count; ++i) {
        h = hash_combine(h, geom::hash(objs[i]));
    }
    return h;
}

inline size_t hash_bytes(const void* obj, size_t size) {
    const char* bytes = reinterpret_cast<const char*>(obj);
    std::string_view s{bytes, size};
    return std::hash<std::string_view>{}(s);
}


} // namespace geom
