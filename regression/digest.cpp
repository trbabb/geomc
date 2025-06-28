#define TEST_MODULE_NAME Digest

// #include <iostream>

#include <random>
#include <iomanip>
#include <gtest/gtest.h>
#include <geomc/Hash.h>

using namespace geom;
using namespace std;

typedef std::mt19937_64 rng_t;

using uint128_t = __uint128_t;

inline std::ostream& operator<<(std::ostream& os, uint128_t value) {
    uint64_t hi = value >> 64;
    uint64_t lo = value & ~uint64_t(0);
    os  << "0x"
        << std::hex
        << std::setw(16) << std::setfill('0') << hi
        << std::setw(16) << std::setfill('0') << lo
        << std::dec;
    return os;
}

template <typename T>
size_t h(const T& v) {
    return std::hash<T>{}(v);
}

template <typename T>
uint128_t d128(const T& v) {
    return geom::hash<T, uint128_t>(v);
}

template <typename T>
size_t dsz(const T& v) {
    return geom::hash<T, size_t>(v);
}

TEST(TEST_MODULE_NAME, verify_string) {
    std::string s0 = "abc";
    std::string_view s1 {s0};
    char s2[4] = "abc";
    const char* s3 = "abc";
    std::string_view s4 {s2, 3};
    std::string_view s5 {s2, 4};
    
    size_t    h_sz  = h(s0);
    uint128_t h_128 = d128(s0);
    
    EXPECT_EQ(h_sz, dsz(s0));
    EXPECT_EQ(h_sz, dsz(s1));
    EXPECT_EQ(h_sz, dsz(s2));
    EXPECT_EQ(h_sz, dsz(s3));
    EXPECT_EQ(h_sz, dsz(s4));
    EXPECT_NE(h_sz, dsz(s5)); // explicitly includes the null terminator!
    
    EXPECT_EQ(h_128, d128(s1));
    EXPECT_EQ(h_128, d128(s2));
    EXPECT_EQ(h_128, d128(s3));
    EXPECT_EQ(h_128, d128(s4));
    EXPECT_NE(h_128, d128(s5)); // explicitly includes the null terminator!
    
    // std::cout << "s0 (string):                " << digest(s0) << " len " << s0.size() << std::endl;
    // std::cout << "s1 (string_view):           " << digest(s1) << " len " << s1.size() << std::endl;
    // std::cout << "s2 (char[4]):               " << digest(s2) << " len " << 4 << std::endl;
    // std::cout << "s3 (const char*):           " << digest(s3) << " len " << 4 << std::endl;
    // std::cout << "s4 (string_view no null):   " << digest(s4) << " len " << s4.size() << std::endl;
    // std::cout << "s5 (string view with null): " << digest(s5) << " len " << s5.size() << std::endl;
}
