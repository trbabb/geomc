/*
 * Xoshiro.h
 *
 */

#include <atomic>
#include <chrono>
#include <geomc/random/Random.h>

namespace geom {


inline uint64_t rot64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}


uint64_t splitmix64(uint64_t* state) {
    uint64_t x = (*state += 0x9E3779B97f4A7C15);
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EB;
    return x ^ (x >> 31);
}


class XoshiroRand: public geom::Random {
public:
    
    XoshiroRand() {
        static std::atomic<uint64_t> uniq = 0;
        uint64_t x = std::chrono::steady_clock::now().time_since_epoch().count();
        x ^= uniq++;
        rseed(x);
    }
    
    XoshiroRand(uint64_t seed) {
        rseed(seed);
    };
    
    inline void rseed(uint64_t seed) {
        uint64_t x = seed;
        s[0] = splitmix64(&x);
        s[1] = splitmix64(&x);
        s[2] = splitmix64(&x);
        s[3] = splitmix64(&x);
    }
    
    inline uint64_t rand64() {
        uint64_t out = rot64(s[1] * 5, 7) * 9;
        uint64_t   t = s[1] << 17;
        
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        
        s[2] ^= t;
        s[3] = rot64(s[3], 45);
        
        return out;
    }
    
    uint32_t rand32() {
        return rand64() >> 32;
    }
    
private:
    
    uint64_t s[4];
};

} // namespace geom
