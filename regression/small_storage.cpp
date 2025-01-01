#define TEST_MODULE_NAME SmallStorage

#include <iostream>
#include <pcg_random.hpp>
#include <gtest/gtest.h>
#include <geomc/SmallStorage.h>

// todo: moar tests
//   - we don't actually verify that element ctors/dtors are run properly

using namespace geom;
using namespace std;

typedef pcg64 rng_t;

using SStor = SmallStorage<int64_t,4>;


void check_okay(const SStor& s, int start, int increment, size_t sz) {
    EXPECT_EQ(s.size(), sz);
    EXPECT_GE(s.capacity(), sz);
    for (size_t i = 0; i < sz; ++i) {
        EXPECT_EQ(s[i], start + i * increment);
    }
}


TEST(TEST_MODULE_NAME, verify_default_ctor) {
    SStor s{};
    EXPECT_EQ(s.size(), 0);
    EXPECT_EQ(s.static_capacity(), s.capacity());
    EXPECT_EQ(s.begin(), s.end());
    
}

TEST(TEST_MODULE_NAME, verify_initializer_list) {
    SStor ss {{1,2,3,4}};
    EXPECT_EQ(ss[0], 1);
    EXPECT_EQ(ss[1], 2);
    EXPECT_EQ(ss[2], 3);
    EXPECT_EQ(ss[3], 4);
    EXPECT_EQ(ss.size(), 4);
    EXPECT_EQ(ss.capacity(), 4);
    
    SStor s1 {{1,2,3,4,5,6,7,8}};
    for (int i = 0; i < 8; ++i) {
        EXPECT_EQ(s1[i], i + 1);
    }
    EXPECT_EQ(s1.size(), 8);
    EXPECT_GE(s1.capacity(), 8);
}

TEST(TEST_MODULE_NAME, verify_copy_assign) {
    SStor ss {{3,2,1}};
    {
        SStor s0 {{0,1,2,3,4,5}};
        ss = s0;
    }
    EXPECT_EQ(ss.size(), 6);
    EXPECT_GE(ss.capacity(), 6);
    for (int i = 0; i < 6; ++i) {
        EXPECT_EQ(ss[i],i);
    }
    
    {
        SStor s0 {{5,4,3}};
        ss = s0;
    }
    EXPECT_EQ(ss.size(), 3);
    EXPECT_GE(ss.capacity(), 3);
    EXPECT_EQ(ss[0], 5);
    EXPECT_EQ(ss[1], 4);
    EXPECT_EQ(ss[2], 3);
}

TEST(TEST_MODULE_NAME, verify_copy_ctor) {
    SStor ssmol {{5,4,3}};
    SStor sbig  {{0,1,2,3,4,5,6,7}};
    {
        SStor s {ssmol};
        EXPECT_EQ(s.size(), ssmol.size());
        EXPECT_EQ(s.capacity(), ssmol.capacity());
        EXPECT_EQ(s[0], 5);
        EXPECT_EQ(s[1], 4);
        EXPECT_EQ(s[2], 3);
    }
    
    {
        SStor s {sbig};
        EXPECT_EQ(s.size(), sbig.size());
        EXPECT_GE(s.capacity(), s.size());
        for (int i = 0; i < s.size(); ++i) {
            EXPECT_EQ(s[i], i);
        }
    }
}

TEST(TEST_MODULE_NAME, verify_move_assign) {
    SStor s;
    {
        SStor s0 {{2,4,6,8,10,12}};
        s = std::move(s0);
        EXPECT_NE(s0.begin(), s.begin());
    }
    check_okay(s, 2, 2, 6);
    {
        SStor s0 {{5,4,3,2}};
        s = std::move(s0);
        EXPECT_NE(s0.begin(), s.begin());
    }
    check_okay(s, 5, -1, 4);
}

TEST(TEST_MODULE_NAME, verify_move_ctor) {
    SStor ssmol {{0,1,2,3}};
    SStor sbig  {{2,4,6,8,10, 12,14,16,18,20}};
    {
        SStor s {std::move(ssmol)};
        EXPECT_NE(s.begin(), ssmol.begin());
        check_okay(s, 0, 1, 4);
    }
    {
        SStor s {std::move(sbig)};
        EXPECT_NE(s.begin(), sbig.begin());
        check_okay(s, 2, 2, 10);
    }
}

TEST(TEST_MODULE_NAME, verify_reserve) {
    SStor s{{0,1,2,3}};
    s.reserve(20);
    check_okay(s, 0, 1, 4);
    EXPECT_GE(s.capacity(), 20);
}

TEST(TEST_MODULE_NAME, verify_push_back) {
    SStor s;
    for (size_t i = 0; i < 33; ++i) {
        s.push_back(i * 2);
        check_okay(s, 0, 2, i + 1);
    }
}
