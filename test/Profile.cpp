/*
 * Profile.cpp
 *
 *  Created on: Apr 16, 2011
 *      Author: tbabb
 */

#include <time.h>
#include <string>
#include <iostream>

#include <boost/iterator/counting_iterator.hpp>

#include "linalg/Vec.h"
#include "linalg/Matrix.h"
#include "function/Path.h"
#include "random/RandomTools.h"
#include "RandomBattery.h"
#include "random/MTRand.h"
#include "random/LCRand.h"
#include "function/PerlinNoise.h"
#include "function/Raster.h"
#include "shape/BinLatticePartition.h"
#include "shape/Trace.h"

//#include "linalg/mtxdetail/SimpleMatrix.h"

#define NUM_PROFILE_CASES 1024;

using namespace geom;
using namespace std;

template <index_t N> 
void fill_unit_vec_array(typename PointType<double,N>::point_t *dst, index_t n){
    RandomVectors<double> rntools = RandomVectors<double>();
    
    for (index_t i = 0; i < n; i++){
        dst[i] = rntools.template unit<N>();
    }
}

template <index_t N> void fill_unit_raw_vec_array(double dst[][N], index_t n){
    RandomVectors<double> rntools = RandomVectors<double>();
    Vec<double,N> v;
    
    for (index_t i = 0; i < n; i++){
        v = rntools.template unit<N>();
        for (index_t axis = 0; axis < N; axis++){
            dst[i][axis] = v[axis];
        }
    }
}

template <index_t N> void fill_range_vec_array(Vec<double, N> *dst, index_t n, Vec<double,N> lo, Vec<double,N> hi){
    RandomVectors<double> rntools = RandomVectors<double>();
    
    for (index_t i = 0; i < n; i++){
        dst[i] = rntools.box(lo, hi);
    }
}

double profile_vec_cross(index_t iters){
    const index_t n = NUM_PROFILE_CASES;
    Vec3d v;
    Vec3d vecs_src[n];
    Vec3d vecs_dst[n];
    
    fill_unit_vec_array<3>(vecs_src, n);
    
    index_t idx = 0;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++){
        vecs_dst[idx % n] = vecs_src[idx % n].cross(vecs_src[(idx+1) % n]);
        idx += 1;
    }
    clock_t end = clock();
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <index_t N> double profile_vec_add(index_t iters){
    index_t n = NUM_PROFILE_CASES;
    Vec<double, N> vecs_src_1[n];
    Vec<double, N> vecs_src_2[n];
    Vec<double, N> vecs_dst[n];
    
    fill_unit_vec_array<N>(vecs_src_1, n);
    fill_unit_vec_array<N>(vecs_src_2, n);
    
    index_t idx = 0;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++){
        vecs_dst[idx] = vecs_src_1[idx] + (vecs_src_2[idx]);
        idx = (idx + 1) % n;
    }
    clock_t end = clock();
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <index_t N> double profile_vec_dot(index_t iters){
    index_t n = NUM_PROFILE_CASES;
    Vec<double, N> vecs_src_1[n];
    Vec<double, N> vecs_src_2[n];
    double vecs_dst[n];
    
    fill_unit_vec_array<N>(vecs_src_1, n);
    fill_unit_vec_array<N>(vecs_src_2, n);
    
    index_t idx = 0;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++){
        vecs_dst[idx] = vecs_src_1[idx].dot(vecs_src_2[idx]);
        idx = (idx + 1) % n;
    }
    clock_t end = clock();
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <index_t N> double profile_vec_norm(index_t iters){
    index_t n = NUM_PROFILE_CASES;
    Vec<double, N> vecs_src[n];
    double vecs_dst[n];
    
    fill_unit_vec_array<N>(vecs_src, n);
    
    index_t idx = 0;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++){
        vecs_dst[idx % n] = vecs_src[idx % n].mag();
        idx += 1;
    }
    clock_t end = clock();
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <index_t N> double profile_vec_norm2(index_t iters){
    index_t n = NUM_PROFILE_CASES;
    Vec<double, N> vecs_src[n];
    double vecs_dst[n];
    
    fill_unit_vec_array<N>(vecs_src, n);
    
    index_t idx = 0;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++){
        vecs_dst[idx % n] = vecs_src[idx % n].mag2();
        idx += 1;
    }
    clock_t end = clock();
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <index_t N> double profile_vec_hash(index_t iters){
    index_t n = NUM_PROFILE_CASES;
    Vec<double, N> vecs_src[n];
    index_t vecs_dst[n];
    
    fill_unit_vec_array<N>(vecs_src, n);
    
    index_t idx = 0;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++){
        vecs_dst[idx % n] = vecs_src[idx % n].hashcode();
        idx += 1;
    }
    clock_t end = clock();
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <index_t N> double profile_raw_vec_add(index_t iters){
    index_t n = NUM_PROFILE_CASES;
    double vecs_src1[n][N];
    double vecs_src2[n][N];
    double vecs_dst[n][N];
    
    fill_unit_raw_vec_array(vecs_src1, n);
    fill_unit_raw_vec_array(vecs_src2, n);
    
    index_t idx = 0;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++){
        for (index_t axis = 0; axis < N; axis++){
            vecs_dst[idx][axis] = vecs_src1[idx][axis] + vecs_src2[idx][axis];
        }
        idx = (idx + 1) % n;
    }
    clock_t end = clock();
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <index_t N> double profile_perlin(index_t iters){
	typedef typename PointType<double,N>::point_t point_t;
    index_t n = NUM_PROFILE_CASES;
    point_t range = point_t(100000);
    point_t vecs_src[n];
    PerlinNoise<double,N> perlin;
    double dest_vals[n];
    
    fill_range_vec_array<N>(vecs_src, n, -range, range);
    
    index_t idx = 0;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++){
        dest_vals[idx] = perlin.eval(vecs_src[n]);
        idx = (idx + 1) % n;
    }
    clock_t end = clock();
    return (end-start) / (double)CLOCKS_PER_SEC;
}

// careful with this. allocates lots o memory.
template <typename T> double profile_rayTriangleTest(index_t iters){
    RandomVectors<T> rvs;
	Ray<T,3> *rays = new Ray<T,3>[iters];
	// generate N random rays. 
	for (index_t i = 0; i < iters; i++){
		rays[i].origin    = rvs.template unit<3>(2.0);
		rays[i].direction = rvs.template unit<3>();
	}
	// triangle points
	Vec<T,3> p0(0.75,0,0);
	Vec<T,3> p1(0,0.75,0);
	Vec<T,3> p2(0,0,0.75);

	clock_t start = clock();
	for (index_t i = 0; i < iters; i++){
		trace_tri(p0, p1, p2, rays[i], HIT_FRONT);
	}
	clock_t end = clock();
	
	delete [] rays;

	return (end-start) / (double)CLOCKS_PER_SEC;
}

template <typename T, index_t N> double profile_path(index_t iters){
    Path<T,N> p;
    RandomVectors<T> rvs;
    int n_knots = 50;
    for (index_t i = 0; i < n_knots; i++){
        p.knots.push_back(Ray<T,N>(rvs.template box<N>(), rvs.template solidball<N>()));
    }
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++){
        p.eval(getRandom()->rand(0,n_knots));
    }
    clock_t end = clock();
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <typename T, index_t N, index_t Channels, Interpolation Interp> 
double profile_raster(index_t iters){
	index_t n = 128;
	Vec<index_t,N> dim(n);
    Raster<T,T,N,Channels> image(dim);
    Vec<T,N> coords[n];
    RandomVectors<T> rvs;
    
    for (index_t i = 0; i < n; i++){
    	coords[i] = rvs.box(Vec<T,N>::zeros, (Vec<T,N>)dim);
    	for (index_t j = 0; j < n; j++){
    		image.set(Vec<index_t,N>(i,j), rvs.template unit<Channels>());
    	}
    }
    
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++){
        image.template sample<EDGE_CLAMP,Interp>(coords[i & (n-1)]);
    }
    clock_t end = clock();
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <typename T> double profile_rotCtr(index_t iters){
    RandomVectors<T> rvs;
    AffineTransform<T,3> xf;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++){
        T angle = rvs.rng->template rand<T>(0, 2*M_PI);
        Vec<T,3> axis = rvs.template unit<3>();
        Vec<T,3> ctr  = rvs.box(Vec<T,3>(-10,-10,-10), Vec<T,3>(10,10,10));
        xf = rotation<T>(axis, ctr, angle);
    }
    clock_t end = clock();
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <typename M> double profile_mtxCopy(index_t iters) {
	index_t rows = 5;
	index_t cols = 5;
	M m(rows,cols);
	Random *rng = getRandom();
	typedef typename M::elem_t T;
	T *ary = new T[rows*cols+1];
	for (index_t i = 0; i < rows*cols+1; i++){
		ary[i] = rng->rand(1.0);
	}
	clock_t start = clock();
	for (index_t i = 0; i < iters; i++){
		std::copy(m.begin(), m.end(), ary);	
	}
	clock_t end = clock();
	
	delete [] ary;
	return (end-start) / (double)CLOCKS_PER_SEC;
}

template <typename M> double profile_mtxRegionCopy(index_t iters) {
	index_t rows = 5;
	index_t cols = 5;
	M m(rows,cols);
	Random *rng = getRandom();
	typedef typename M::elem_t T;
	T *ary = new T[rows*cols+1];
	for (index_t i = 0; i < rows*cols+1; i++){
		ary[i] = rng->rand(1.0);
	}
	clock_t start = clock();
	for (index_t i = 0; i < iters; i++){
		MatrixRegion region(MatrixCoord(1,1), MatrixCoord(rows-1, cols-1));
		std::copy(m.region_begin(region), m.region_end(region), ary);
	}
	clock_t end = clock();
	
	delete [] ary;
	return (end-start) / (double)CLOCKS_PER_SEC;
}

template <typename T, index_t N> double profile_mtxInverse(index_t iters) {
	SimpleMatrix<T,N,N> mtx[2];
	Random *rng = getRandom();
	for (T *p = mtx[0].begin(); p != mtx[0].end(); p++){
		*p = rng->rand(1.0);
	}
	index_t x = 0;
	clock_t start = clock();
	for (index_t i = 0; i < iters; i++, x = !x){
		inv(&mtx[!x], mtx[x]);	
	}
	clock_t end = clock();
	
	return (end-start) / (double)CLOCKS_PER_SEC;
}

template <typename T, index_t N> double profile_mtxInverseLU(index_t iters) {
	SimpleMatrix<T,N,N> mtx[2];
	Random *rng = getRandom();
	for (T *p = mtx[0].begin(); p != mtx[0].end(); p++){
		*p = rng->rand(1.0);
	}
	index_t x = 0;
	clock_t start = clock();
	for (index_t i = 0; i < iters; i++, x = !x){
		PLUDecomposition<T,N,N> plu = plu_decompose(mtx[x]);
		plu.inverse(&mtx[!x]);	
	}
	clock_t end = clock();
	
	return (end-start) / (double)CLOCKS_PER_SEC;
}

template <typename T, index_t N> void test_mtxInverse(index_t iters) {
	SimpleMatrix<T,N,N> mtx[2];
	Random *rng = getRandom();
	for (T *p = mtx[0].begin(); p != mtx[0].end(); p++){
		*p = rng->rand(1.0);
	}
	std::cout << mtx[0];
	index_t x = 0;
	for (index_t i = 0; i < iters; i++, x = !x){
		inv(&mtx[!x], mtx[x]);
		std::cout << mtx[!x];
	}
}

template <typename M> void test_mtxRegionCopy(index_t rows, index_t cols) {
	M m(rows,cols);
	typedef typename M::elem_t T;
	const index_t subn = (rows-2)*(cols-2);
	T *v = new T[subn];
	
	for (index_t i = 0; i < subn; i++){
		v[i] = 11*i;
	}
	
	std::fill(m.begin(), m.end(), 1);
	std::copy(v, v+subn, m.region_begin(1,1,rows-1, cols-1));
	cout << m;
}

template <index_t N> void test_permuteMatrix() {
	PermutationMatrix<N> p;
	Vec<double,N> v;
	Random *rng = getRandom();
	for (index_t i = 0; i < N; i++) {
		v[i] = i;
	}
	for (index_t i = 0; i < N*2; i++) {
		p.swap_rows(rng->rand(N), rng->rand(N));
	}
	SimpleMatrix<double,N,N> m;
	SimpleMatrix<double,N,N> minv;
	std::copy(p.begin(), p.end(), m.begin());
	inv(&minv, p);
	cout << "v: " << v << endl;
	cout << "P: " << endl << p << endl;
	cout << "P*v: " << (p * v) << endl;
	cout << "M: " << endl << m << endl;
	cout << "M*v: " << (m * v) << endl;
	cout << "P^-1: " << endl << minv << endl;
}

template <index_t M, index_t N> void test_simpleMatrix() {
	SimpleMatrix<double, M, N> m;
	SimpleMatrix<double, DYNAMIC_DIM, DYNAMIC_DIM> q(M,N);
	m[0][1] = 5;
	q[0][1] = 8;
	cout << m << endl;
	cout << q << endl;
}

void test_matrixOrder() {
	index_t sz_r = 5;
	index_t sz_c = 7;
	SparseMatrix<double> mtx(sz_r,sz_c);

	double ct = 0;
	for (SparseMatrix<double>::iterator it = mtx.begin(); it != mtx.end(); it++) {
		*it = ct++;
	}
    
	std::cout << "[";
	for (index_t r = 0; r < sz_r; r++) {
		for (index_t c = 0; c < sz_c; c++) {
			std::cout << mtx.get(r,c) << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "]" << std::endl;
}

void profile(std::string name, double (*fnc)(index_t), index_t iterations){
    double t = fnc(iterations);
    std::cout << "profile " << name << " (" << iterations << " iters): " << t << " secs; " << iterations/t << " ops/sec"  << std::endl;
}

void test_matrix_asplode(){
    SimpleMatrix<double,4,4> mtx;
    SimpleMatrix<double,2,5> other;
    // correctly, won't compile:
    // mtx * other;
}

using geom::operator<<;

//int main_profile(int arc, char **argv){
int main(int argc, char** argv) {
    Vec2d a2d;
    Vec3d a3d;
    Vec4d a4d;
    
    index_t iters = 10000000;
    
    test_matrixOrder();
    
    assert(&a2d.get(0) == &a2d.x && &a2d.get(1) == &a2d.y);
    assert(&a3d.get(0) == &a3d.x && &a3d.get(1) == &a3d.y && &a3d.get(2) == &a3d.z);
    assert(&a4d.get(0) == &a4d.x && &a4d.get(1) == &a4d.y && &a4d.get(2) == &a4d.z && &a4d.get(3) == &a4d.w);
    
    test_permuteMatrix<5>();
    test_simpleMatrix<4,3>();
    
    profile("3d cross product", profile_vec_cross, iters);
    std::cout << std::endl;
    
    profile("2d add", profile_vec_add<2>, iters);
    profile("3d add", profile_vec_add<3>, iters);
    profile("4d add", profile_vec_add<4>, iters);
    profile("8d add", profile_vec_add<8>, iters);
    std::cout << std::endl;
    
    profile("2d raw add", profile_raw_vec_add<2>, iters);
    profile("3d raw add", profile_raw_vec_add<3>, iters);
    profile("4d raw add", profile_raw_vec_add<4>, iters);
    profile("8d raw add", profile_raw_vec_add<8>, iters);
    std::cout << std::endl;
    
    profile("2d norm", profile_vec_norm<2>, iters);
    profile("3d norm", profile_vec_norm<3>, iters);
    profile("4d norm", profile_vec_norm<4>, iters);
    profile("8d norm", profile_vec_norm<8>, iters);
    std::cout << std::endl;
    
    profile("2d norm2", profile_vec_norm2<2>, iters);
    profile("3d norm2", profile_vec_norm2<3>, iters);
    profile("4d norm2", profile_vec_norm2<4>, iters);
    profile("8d norm2", profile_vec_norm2<8>, iters);
    std::cout << std::endl;
    
    //profile("1d perlin", profile_perlin<1>, iters/100);
    profile("2d perlin", profile_perlin<2>, iters/100);
    profile("3d perlin", profile_perlin<3>, iters/100);
    profile("4d perlin", profile_perlin<4>, iters/100);
    profile("5d perlin", profile_perlin<5>, iters/100);
    std::cout << std::endl;
    
    profile("2d hash", profile_vec_hash<2>, iters);
    profile("3d hash", profile_vec_hash<3>, iters);
    profile("4d hash", profile_vec_hash<4>, iters);
    std::cout << std::endl;
    
    profile("2d dot", profile_vec_dot<2>, iters);
    profile("3d dot", profile_vec_dot<3>, iters);
    profile("4d dot", profile_vec_dot<4>, iters);
    profile("8d dot", profile_vec_dot<8>, iters);
    std::cout << std::endl;
	
    profile("rayf-triangle hit", profile_rayTriangleTest<float>,  iters/10);
    profile("rayd-triangle hit", profile_rayTriangleTest<double>, iters/10);
	std::cout << std::endl;
    
    profile("2d path", profile_path<double, 2>, iters);
    profile("3d path", profile_path<double, 3>, iters);
    profile("4d path", profile_path<double, 4>, iters);
    profile("2f path", profile_path<float, 2>, iters);
    profile("3f path", profile_path<float, 3>, iters);
    profile("4f path", profile_path<float, 4>, iters);
    std::cout << std::endl;
    
    profile("2f->3f linear img sample", profile_raster<float, 2, 3, INTERP_LINEAR>, iters/100);
    profile("3f->3f linear img sample", profile_raster<float, 3, 3, INTERP_LINEAR>, iters/100);
    profile("2f->3f cubic img sample",  profile_raster<float, 2, 3, INTERP_CUBIC>,  iters/100);
    profile("3f->3f cubic img sample",  profile_raster<float, 3, 3, INTERP_CUBIC>,  iters/100);
    std::cout << std::endl;
    
    profile("5x5 mtxf region copy", profile_mtxRegionCopy<SimpleMatrix<float,0,0> >,    iters);
    profile("5x5 mtxd region copy", profile_mtxRegionCopy<SimpleMatrix<double,0,0> >,   iters);
    profile("5x5 sparsef region copy", profile_mtxRegionCopy<SparseMatrix<float> >,  iters);
    profile("5x5 sparsed region copy", profile_mtxRegionCopy<SparseMatrix<double> >, iters);
    profile("5x5 diagf region copy", profile_mtxRegionCopy<DiagMatrix<float,5,5> >, iters);
    std::cout << std::endl;

    profile("5x5 mtxf copy", profile_mtxCopy<SimpleMatrix<float,0,0> >, iters/100);
    profile("5x5 mtxd copy", profile_mtxCopy<SimpleMatrix<double,0,0> >, iters/100);
    profile("5x5 sparsef copy", profile_mtxCopy<SparseMatrix<float> >, iters/100);
    profile("5x5 sparsed copy", profile_mtxCopy<SparseMatrix<double> >, iters/100);
    profile("5x5 diagf copy", profile_mtxCopy<DiagMatrix<float,5,5> >, iters/100);
    std::cout << std::endl;
    
    profile("2x2 mtxf inv", profile_mtxInverse<float,2>,  iters/10);
    profile("3x3 mtxf inv", profile_mtxInverse<float,3>,  iters/10);
    profile("4x4 mtxf inv", profile_mtxInverse<float,4>,  iters/10);
    profile("2x2 mtxd inv", profile_mtxInverse<double,2>, iters/10);
    profile("3x3 mtxd inv", profile_mtxInverse<double,3>, iters/10);
    profile("4x4 mtxd inv", profile_mtxInverse<double,4>, iters/10);
    profile("8x8 mtxd inv", profile_mtxInverse<double,8>, iters/100);
    std::cout << std::endl;
    
    profile("rotationd matrix", profile_rotCtr<double>, iters/10);
    profile("rotationf matrix", profile_rotCtr<float>, iters/10);
    std::cout << std::endl;

    std::cout << "sizeof vec2f: " << sizeof(Vec<float,2>) << std::endl;
    std::cout << "sizeof vec3f: " << sizeof(Vec<float,3>) << std::endl;
    std::cout << "sizeof vec4f: " << sizeof(Vec<float,4>) << std::endl;
    std::cout << "sizeof vec2d: " << sizeof(Vec<double,2>) << std::endl;
    std::cout << "sizeof vec3d: " << sizeof(Vec<double,3>) << std::endl;
    std::cout << "sizeof vec4d: " << sizeof(Vec<double,4>) << std::endl;
    
    try {
        test_matrix_asplode();
    } catch (GeomException &ex){
        std::cout << ex.what() << std::endl << std::endl;
    }
    
    /*
    std::cout << "Mersenne Twister" << std::endl << "================================" << std::endl;
    RandomBattery bat = RandomBattery(new MTRand());
    bat.runAll(iters);
    std::cout << std::endl;
    
    std::cout << "Linear Congruential" << std::endl << "================================" << std::endl;
    bat = RandomBattery(new LCRand());
    bat.runAll(iters);
    std::cout << std::endl;
    */
    
    std::cout << "done." << std::endl;
    
    return 0;
}
