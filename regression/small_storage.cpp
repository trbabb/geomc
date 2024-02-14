#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SmallStorage

#include <iostream>
#include <random>
#include <pcg_random.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <geomc/SmallStorage.h>

// todo: moar tests
//   - we don't actually verify that element ctors/dtors are run properly

using namespace geom;
using namespace std;

typedef pcg64 rng_t;

using SStor = SmallStorage<int64_t,4>;


void check_okay(const SStor& s, int start, int increment, size_t sz) {
    BOOST_CHECK_EQUAL(s.size(), sz);
    BOOST_CHECK_GE(s.capacity(), sz);
    for (size_t i = 0; i < sz; ++i) {
        BOOST_CHECK_EQUAL(s[i], start + i * increment);
    }
}


BOOST_AUTO_TEST_CASE(verify_default_ctor) {
    SStor s{};
    BOOST_CHECK_EQUAL(s.size(), 0);
    BOOST_CHECK_EQUAL(s.static_capacity(), s.capacity());
    BOOST_CHECK_EQUAL(s.begin(), s.end());
    
}

BOOST_AUTO_TEST_CASE(verify_initializer_list) {
    SStor ss {{1,2,3,4}};
    BOOST_CHECK_EQUAL(ss[0], 1);
    BOOST_CHECK_EQUAL(ss[1], 2);
    BOOST_CHECK_EQUAL(ss[2], 3);
    BOOST_CHECK_EQUAL(ss[3], 4);
    BOOST_CHECK_EQUAL(ss.size(), 4);
    BOOST_CHECK_EQUAL(ss.capacity(), 4);
    
    SStor s1 {{1,2,3,4,5,6,7,8}};
    for (int i = 0; i < 8; ++i) {
        BOOST_CHECK_EQUAL(s1[i], i + 1);
    }
    BOOST_CHECK_EQUAL(s1.size(), 8);
    BOOST_CHECK_GE(s1.capacity(), 8);
}

BOOST_AUTO_TEST_CASE(verify_copy_assign) {
    SStor ss {{3,2,1}};
    {
        SStor s0 {{0,1,2,3,4,5}};
        ss = s0;
    }
    BOOST_CHECK_EQUAL(ss.size(), 6);
    BOOST_CHECK_GE(ss.capacity(), 6);
    for (int i = 0; i < 6; ++i) {
        BOOST_CHECK_EQUAL(ss[i],i);
    }
    
    {
        SStor s0 {{5,4,3}};
        ss = s0;
    }
    BOOST_CHECK_EQUAL(ss.size(), 3);
    BOOST_CHECK_GE(ss.capacity(), 3);
    BOOST_CHECK_EQUAL(ss[0], 5);
    BOOST_CHECK_EQUAL(ss[1], 4);
    BOOST_CHECK_EQUAL(ss[2], 3);
}

BOOST_AUTO_TEST_CASE(verify_copy_ctor) {
    SStor ssmol {{5,4,3}};
    SStor sbig  {{0,1,2,3,4,5,6,7}};
    {
        SStor s {ssmol};
        BOOST_CHECK_EQUAL(s.size(), ssmol.size());
        BOOST_CHECK_EQUAL(s.capacity(), ssmol.capacity());
        BOOST_CHECK_EQUAL(s[0], 5);
        BOOST_CHECK_EQUAL(s[1], 4);
        BOOST_CHECK_EQUAL(s[2], 3);
    }
    
    {
        SStor s {sbig};
        BOOST_CHECK_EQUAL(s.size(), sbig.size());
        BOOST_CHECK_GE(s.capacity(), s.size());
        for (int i = 0; i < s.size(); ++i) {
            BOOST_CHECK_EQUAL(s[i], i);
        }
    }
}

BOOST_AUTO_TEST_CASE(verify_move_assign) {
    SStor s;
    {
        SStor s0 {{2,4,6,8,10,12}};
        s = std::move(s0);
        BOOST_CHECK_NE(s0.begin(), s.begin());
    }
    check_okay(s, 2, 2, 6);
    {
        SStor s0 {{5,4,3,2}};
        s = std::move(s0);
        BOOST_CHECK_NE(s0.begin(), s.begin());
    }
    check_okay(s, 5, -1, 4);
}

BOOST_AUTO_TEST_CASE(verify_move_ctor) {
    SStor ssmol {{0,1,2,3}};
    SStor sbig  {{2,4,6,8,10, 12,14,16,18,20}};
    {
        SStor s {std::move(ssmol)};
        BOOST_CHECK_NE(s.begin(), ssmol.begin());
        check_okay(s, 0, 1, 4);
    }
    {
        SStor s {std::move(sbig)};
        BOOST_CHECK_NE(s.begin(), sbig.begin());
        check_okay(s, 2, 2, 10);
    }
}

BOOST_AUTO_TEST_CASE(verify_reserve) {
    SStor s{{0,1,2,3}};
    s.reserve(20);
    check_okay(s, 0, 1, 4);
    BOOST_CHECK_GE(s.capacity(), 20);
}

BOOST_AUTO_TEST_CASE(verify_push_back) {
    SStor s;
    for (size_t i = 0; i < 33; ++i) {
        s.push_back(i * 2);
        check_okay(s, 0, 2, i + 1);
    }
}
