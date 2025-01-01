// void destructive_apply_permutation(index_t* p, T* A, index_t n, index_t m);

#define TEST_MODULE_NAME LuDebug

#include <gtest/gtest.h>
#include <geomc/linalg/Matrix.h>
#include <geomc/random/RandomTools.h>


using namespace geom;
using namespace std;


TEST(TEST_MODULE_NAME, verify_permute) {
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
    EXPECT_EQ(m0, m1);
}
