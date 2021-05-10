#include <chrono>
#include <algorithm>

#include "shapesample.h"

using namespace geom;
using namespace std;
using namespace std::chrono;


template <typename Shape>
struct PerfTest {
    static void test(rng_t* rng, index_t ct) {}
};


template <typename T, index_t N>
struct PerfTest <Simplex<T,N>> {
    static T test(rng_t* rng, Simplex<T,N>* shapes, index_t ct) {
        T x = 0; // prevent optimizing out the loop body
        for (index_t i = 0; i < ct; ++i) {
            auto& s   = shapes[i];
            auto bnd  = s.bounds();
            auto ctr  = bnd.center();
            auto dims = bnd.dimensions();
            for (index_t j = 0; j < 100; ++j) {
                auto p = rnd<T,N>(rng) * dims + ctr;
                x += s.project(p)[0];
            }
        }
        return x;
    }
};


template <typename Shape>
void fill_shape_buffer(rng_t* rng, Shape* shapes, index_t ct) {
    for (index_t i = 0; i < ct; ++i) {
        shapes[i] = RandomShape<Shape>::rnd_shape(rng);
    }
}


template <typename Shape>
void time_perftest_shape(index_t ct) {
    rng_t rng(18374691138699945602ULL);
    index_t batch_sz = std::min<index_t>((1<<13), ct);
    Shape* shapes = new Shape[batch_sz];
    typename high_resolution_clock::duration accum_time(0);
    
    for (index_t n = ct; n > 0; n -= batch_sz) {
        index_t cur_batch_sz = std::min(batch_sz, n);
        fill_shape_buffer<Shape>(&rng, shapes, cur_batch_sz);
        auto t0 = high_resolution_clock::now();
        PerfTest<Shape>::test(&rng, shapes, cur_batch_sz);
        auto t1 = high_resolution_clock::now();
        accum_time += t1 - t0;
    }
    
    auto secs  = duration_cast<duration<double>>(accum_time).count();
    std::cout << "perftest " << typeid(Shape).name() << ": ";
    std::cout << (ct / secs) << " ops/s (" << secs << " s)\n";
    delete [] shapes;
}

int main(int argc, char** argv) {
    time_perftest_shape<Simplex<double,2>>(4096);
    time_perftest_shape<Simplex<double,3>>(4096);
    time_perftest_shape<Simplex<double,4>>(4096);
    time_perftest_shape<Simplex<double,5>>(4096);
    time_perftest_shape<Simplex<double,6>>(4096);
    time_perftest_shape<Simplex<double,7>>(4096);
}