#pragma once

#include <atomic>
#include <chrono>
#include <geomc/random/Random.h>

namespace geom {


inline uint64_t splitmix64(uint64_t* state) {
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
    
    XoshiroRand(const uint64_t entropy[4]) {
        seed(entropy);
    }
    
    void rseed(uint64_t seed) {
        uint64_t x = seed;
        s[0] = splitmix64(&x);
        s[1] = splitmix64(&x);
        s[2] = splitmix64(&x);
        s[3] = splitmix64(&x);
    }
    
    template <typename F>
    void seed(F&& entropy) {
        // wangjangle the entropy a bit; this ensures (e.g.) that there are no zeros
        uint64_t x = entropy();
        s[0] = splitmix64(&x);
        x = entropy();
        s[1] = splitmix64(&x);
        x = entropy();
        s[2] = splitmix64(&x);
        x = entropy();
        s[3] = splitmix64(&x);
    }
    
    void seed(const uint64_t entropy[4]) {
        s[0] = entropy[0];
        s[1] = entropy[1];
        s[2] = entropy[2];
        s[3] = entropy[3];
    }
    
    uint64_t rand64() {
        // the xoshiro256++ algorithm:
        uint64_t out = std::rotl<uint64_t>(s[0] + s[3], 23) + s[0];
        uint64_t   t = s[1] << 17;
        
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        
        s[2] ^= t;
        s[3] = std::rotl<uint64_t>(s[3], 45);
        
        return out;
    }
    
    uint32_t rand32() {
        return rand64() >> 32;
    }
    
private:
    
    uint64_t s[4];
};

} // namespace geom
