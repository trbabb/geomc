#include <geomc/geomc_defs.h>

namespace geom {

static inline constexpr void u32to8(uint8_t p[4], uint32_t v) {
    p[0] = (uint8_t)(v);
    p[1] = (uint8_t)(v >> 8);
    p[2] = (uint8_t)(v >> 16);
    p[3] = (uint8_t)(v >> 24);
}

static inline constexpr void u64to8(uint8_t p[8], uint64_t v) {
    u32to8(p, v);
    u32to8(p + 4, v >> 32);
}

static inline constexpr uint64_t u8to64(const uint8_t p[8]) {
    uint64_t v = p[0];
    v |= ((uint64_t) p[1] <<  8) | ((uint64_t) p[2] << 16);
    v |= ((uint64_t) p[3] << 24) | ((uint64_t) p[4] << 32);
    v |= ((uint64_t) p[5] << 40) | ((uint64_t) p[6] << 48);
    v |= ((uint64_t) p[7] << 56);
    return v;
}

static inline void sipround(uint64_t& v0, uint64_t& v1, uint64_t& v2, uint64_t& v3) {
    v0 += v1;
    v1 = std::rotl<uint64_t>(v1, 13);
    v1 ^= v0;
    v0 = std::rotl<uint64_t>(v0, 32);
    v2 += v3;
    v3 = std::rotl<uint64_t>(v3, 16);
    v3 ^= v2;
    v0 += v3;
    v3 = std::rotl<uint64_t>(v3, 21);
    v3 ^= v0;
    v2 += v1;
    v1 = std::rotl<uint64_t>(v1, 17);
    v1 ^= v2;
    v2 = std::rotl<uint64_t>(v2, 32);
}


namespace detail {

/*
    Computes a SipHash value
    *in:    pointer to input data (read-only)
    inlen:  input data length in bytes (any size_t value)
    *k:     pointer to the key data (read-only), must be 16 bytes 
    *out:   pointer to output data (write-only), outlen bytes must be allocated
    outlen: length of the output in bytes, must be 8 or 16
*/
void siphash(
        const void*    in,
        const size_t   inlen,
        const void*    k,
              uint8_t* out,
        const size_t   outlen)
{
    /* default: SipHash-2-4 */
    constexpr size_t c_rounds = 2;
    constexpr size_t d_rounds = 4;
    
    const unsigned char *ni = (const unsigned char *)in;
    const unsigned char *kk = (const unsigned char *)k;
    
    uint64_t v0 = 0x736f6d6570736575;
    uint64_t v1 = 0x646f72616e646f6d;
    uint64_t v2 = 0x6c7967656e657261;
    uint64_t v3 = 0x7465646279746573;
    uint64_t k0 = u8to64(kk);
    uint64_t k1 = u8to64(kk + 8);
    uint64_t m;
    int i;
    const unsigned char *end = ni + inlen - (inlen % sizeof(uint64_t));
    const int left = inlen & 7;
    uint64_t b = ((uint64_t)inlen) << 56;
    v3 ^= k1;
    v2 ^= k0;
    v1 ^= k1;
    v0 ^= k0;
    
    if (outlen == 16)
        v1 ^= 0xee;
    
    for (; ni != end; ni += 8) {
        m = u8to64(ni);
        v3 ^= m;
        for (i = 0; i < c_rounds; ++i) {
            sipround(v0, v1, v2, v3);
        }
        v0 ^= m;
    }
    
    switch (left) {
        case 7: b |= ((uint64_t)ni[6]) << 48; [[fallthrough]];
        case 6: b |= ((uint64_t)ni[5]) << 40; [[fallthrough]];
        case 5: b |= ((uint64_t)ni[4]) << 32; [[fallthrough]];
        case 4: b |= ((uint64_t)ni[3]) << 24; [[fallthrough]];
        case 3: b |= ((uint64_t)ni[2]) << 16; [[fallthrough]];
        case 2: b |= ((uint64_t)ni[1]) <<  8; [[fallthrough]];
        case 1: b |= ((uint64_t)ni[0]); break;
        case 0:
            break;
    }
    
    v3 ^= b;
    for (i = 0; i < c_rounds; ++i) sipround(v0, v1, v2, v3);
    v0 ^= b;
    if (outlen == 16) {
        v2 ^= 0xee;
    } else {
        v2 ^= 0xff;
    }
    for (i = 0; i < d_rounds; ++i) sipround(v0, v1, v2, v3);
    b = v0 ^ v1 ^ v2 ^ v3;
    u64to8(out, b);
    if (outlen == 8) return;
    v1 ^= 0xdd;
    for (i = 0; i < d_rounds; ++i) sipround(v0, v1, v2, v3);
    
    b = v0 ^ v1 ^ v2 ^ v3;
    u64to8(out + 8, b);
    
    return;
}

} // namespace detail
} // namespace geom
