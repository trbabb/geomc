// void destructive_apply_permutation(index_t* p, T* A, index_t n, index_t m);

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE LuDebug

#include <boost/test/unit_test.hpp>
#include <geomc/linalg/Matrix.h>
#include <geomc/random/RandomTools.h>


using namespace geom;
using namespace std;

BOOST_AUTO_TEST_SUITE(lu_debug)

BOOST_AUTO_TEST_CASE(verify_permute) {
    Random* rng = getRandom();
    constexpr index_t N = 5;
    index_t P[N];
    SimpleMatrix<float,N,3> m;
    SimpleMatrix<float,N,3> m0;
    SimpleMatrix<float,N,3> m1;
    PermutationMatrix<N> Px;
    for (index_t i = 0; i < N; ++i) {
        P[i] = i;
        m(i,0) = rng->rand<float>();
        m(i,1) = rng->rand<float>();
        m(i,2) = rng->rand<float>();
    }
    permute(P, N);
    Px.setRowSources(P);
    mul(&m0, Px, m);
    
    mtxcopy(&m1, m);
    detail::destructive_apply_permutation<float,true>(P, m1.data_begin(), N, 3);
    BOOST_CHECK_EQUAL(m0, m1);
}

BOOST_AUTO_TEST_SUITE_END()