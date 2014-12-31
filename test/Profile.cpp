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

#include <geomc/function/Dual.h>
#include <geomc/linalg/Vec.h>
#include <geomc/linalg/Matrix.h>
#include <geomc/random/RandomTools.h>
#include <geomc/random/MTRand.h>
#include <geomc/random/LCRand.h>
#include <geomc/function/Path.h>
#include <geomc/function/PerlinNoise.h>
#include <geomc/function/Raster.h>
#include <geomc/function/SphericalHarmonics.h>
#include <geomc/shape/BinLatticePartition.h>
#include <geomc/shape/Trace.h>
#include <geomc/shape/OrientedRect.h>
#include <geomc/shape/Intersect.h>

#ifdef ENABLE_FRUSTUM
#include <geomc/shape/Frustum.h>
#endif

#include "RandomBattery.h"

#define NUM_PROFILE_CASES 1024;

using namespace geom;
using namespace std;

template <typename T, index_t N> 
void fill_unit_vec_array(typename PointType<T,N>::point_t *dst, index_t n) {
    Sampler<T> rntools = Sampler<T>();
    
    for (index_t i = 0; i < n; i++) {
        dst[i] = rntools.template unit<N>();
    }
}

template <index_t N> void fill_unit_raw_vec_array(double dst[][N], index_t n) {
    Sampler<double> rntools = Sampler<double>();
    Vec<double,N> v;
    
    for (index_t i = 0; i < n; i++) {
        v = rntools.template unit<N>();
        for (index_t axis = 0; axis < N; axis++) {
            dst[i][axis] = v[axis];
        }
    }
}

template <typename T, index_t N> void fill_range_vec_array(Vec<T, N> *dst, index_t n, Vec<T,N> lo, Vec<T,N> hi) {
    Sampler<double> rntools = Sampler<double>();
    
    for (index_t i = 0; i < n; i++) {
        dst[i] = rntools.box(lo, hi);
    }
}

template <typename T, index_t N> void fill_range_vec_array(Vec<Dual<T>, N> *dst, index_t n, Vec<Dual<T>,N> lo, Vec<Dual<T>,N> hi) {
    Sampler< Dual<T> > rntools = Sampler< Dual<T> >();
    
    for (index_t i = 0; i < n; i++) {
        dst[i] = rntools.box(lo, hi);
    }
}

double profile_vec_cross(index_t iters) {
    const index_t n = NUM_PROFILE_CASES;
    Vec3d v;
    Vec3d *vecs_src = new Vec3d[n];
    Vec3d *vecs_dst = new Vec3d[n];
    
    fill_unit_vec_array<double,3>(vecs_src, n);
    
    index_t idx = 0;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++) {
        vecs_dst[idx % n] = vecs_src[idx % n].cross(vecs_src[(idx+1) % n]);
        idx += 1;
    }
    clock_t end = clock();
    
    delete [] vecs_src;
    delete [] vecs_dst;
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <index_t N> double profile_vec_add(index_t iters) {
    index_t n = NUM_PROFILE_CASES;
    Vec<double, N> *vecs_src_1 = new Vec<double,N>[n];
    Vec<double, N> *vecs_src_2 = new Vec<double,N>[n];
    Vec<double, N> *vecs_dst   = new Vec<double,N>[n];
    
    fill_unit_vec_array<double,N>(vecs_src_1, n);
    fill_unit_vec_array<double,N>(vecs_src_2, n);
    
    index_t idx = 0;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++) {
        vecs_dst[idx] = vecs_src_1[idx] + (vecs_src_2[idx]);
        idx = (idx + 1) % n;
    }
    clock_t end = clock();
    
    delete [] vecs_src_1;
    delete [] vecs_src_2;
    delete [] vecs_dst;
    
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <index_t N> double profile_vec_dot(index_t iters) {
    index_t n = NUM_PROFILE_CASES;
    Vec<double, N> *vecs_src_1 = new Vec<double,N>[n];
    Vec<double, N> *vecs_src_2 = new Vec<double,N>[n];
    double *vecs_dst = new double[n];
    
    fill_unit_vec_array<double,N>(vecs_src_1, n);
    fill_unit_vec_array<double,N>(vecs_src_2, n);
    
    index_t idx = 0;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++) {
        vecs_dst[idx] = vecs_src_1[idx].dot(vecs_src_2[idx]);
        idx = (idx + 1) % n;
    }
    clock_t end = clock();
    vecs_dst[0] += 1;
    
    delete [] vecs_src_1;
    delete [] vecs_src_2;
    delete [] vecs_dst;
    
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <index_t N> double profile_vec_norm(index_t iters) {
    index_t n = NUM_PROFILE_CASES;
    Vec<double, N> *vecs_src = new Vec<double,N>[n];
    double *vecs_dst = new double[n];
    
    fill_unit_vec_array<double,N>(vecs_src, n);
    
    index_t idx = 0;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++) {
        vecs_dst[idx % n] = vecs_src[idx % n].mag();
        idx += 1;
    }
    clock_t end = clock();
    vecs_dst[0] += 1;
    
    delete [] vecs_src;
    delete [] vecs_dst;
    
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <index_t N> double profile_vec_norm2(index_t iters) {
    index_t n = NUM_PROFILE_CASES;
    Vec<double, N> *vecs_src = new Vec<double,N>[n];
    double *vecs_dst = new double[n];
    
    fill_unit_vec_array<double,N>(vecs_src, n);
    
    index_t idx = 0;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++) {
        vecs_dst[idx % n] = vecs_src[idx % n].mag2();
        idx += 1;
    }
    clock_t end = clock();
    vecs_dst[0] += 1; //shenanigans to prevent compiler from optimizing to a null loop.
    
    delete [] vecs_src;
    delete [] vecs_dst;
    
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <index_t N> double profile_vec_hash(index_t iters) {
    index_t n = NUM_PROFILE_CASES;
    Vec<double, N> *vecs_src = new Vec<double,N>[n];
    index_t *vecs_dst= new index_t[n];
    
    fill_unit_vec_array<double,N>(vecs_src, n);
    
    index_t idx = 0;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++) {
        vecs_dst[idx % n] = vecs_src[idx % n].hashcode();
        idx += 1;
    }
    clock_t end = clock();
    vecs_dst[0] += 1;
    
    delete [] vecs_src;
    delete [] vecs_dst;
    
    return (end-start) / (double)CLOCKS_PER_SEC;
}

// XXX: TODO: HANDLE THIS

template <index_t N> double profile_raw_vec_add(index_t iters) {
    index_t n = NUM_PROFILE_CASES;
    double vecs_src1[n][N];
    double vecs_src2[n][N];
    double vecs_dst[n][N];
    
    fill_unit_raw_vec_array(vecs_src1, n);
    fill_unit_raw_vec_array(vecs_src2, n);
    
    index_t idx = 0;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++) {
        for (index_t axis = 0; axis < N; axis++) {
            vecs_dst[idx][axis] = vecs_src1[idx][axis] + vecs_src2[idx][axis];
        }
        idx = (idx + 1) % n;
    }
    clock_t end = clock();
    vecs_dst[0][0] += 1;
    
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <index_t N> double profile_perlin(index_t iters) {
	typedef typename PointType<double,N>::point_t point_t;
    index_t n = NUM_PROFILE_CASES;
    point_t range = point_t(100000);
    point_t *vecs_src = new point_t[n];
    PerlinNoise<double,N> perlin;
    double *dest_vals = new double[n];
    
    fill_range_vec_array<double,N>(vecs_src, n, -range, range);
    
    index_t idx = 0;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++) {
        dest_vals[idx] = perlin.eval(vecs_src[i % n]);
        idx = (idx + 1) % n;
    }
    clock_t end = clock();
    dest_vals[0] += 1;
    
    delete [] vecs_src;
    delete [] dest_vals;
    
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <typename T, index_t N> double profile_perlin_grad(index_t iters) {
    typedef Dual<T> dual;
	typedef typename PointType<dual,N>::point_t point_t;
    
    index_t n = NUM_PROFILE_CASES;
    point_t range = point_t(10000);
    point_t *vecs_src = new point_t[n];
    PerlinNoise<dual,N> perlin;
    dual *dest_vals = new dual[n];
    
    fill_range_vec_array<T,N>(vecs_src, n, -range, range);
    
    index_t idx = 0;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++) {
        dest_vals[idx] = perlin.eval(vecs_src[i%n]);
        idx = (idx + 1) % n;
    }
    clock_t end = clock();
    
    delete [] vecs_src;
    delete [] dest_vals;
    
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <typename T, index_t N> double profile_perlin_hard_grad(index_t iters) {
	typedef typename PointType<T,N>::point_t point_t;
    
    index_t n = NUM_PROFILE_CASES;
    point_t range = point_t(10000);
    point_t *vecs_src = new point_t[n];
    PerlinNoise<T,N> perlin;
    point_t *dest_vals = new point_t[n];
    
    fill_range_vec_array<T,N>(vecs_src, n, -range, range);
    
    index_t idx = 0;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++) {
        dest_vals[idx] = perlin.gradient(vecs_src[n]);
        idx = (idx + 1) % n;
    }
    clock_t end = clock();
    
    delete [] vecs_src;
    delete [] dest_vals;
    
    return (end-start) / (double)CLOCKS_PER_SEC;
}

// careful with this. allocates lots o memory.
template <typename T> double profile_rayTriangleTest(index_t iters) {
    Sampler<T> rvs;
	Ray<T,3> *rays = new Ray<T,3>[iters];
	// generate N random rays. 
	for (index_t i = 0; i < iters; i++) {
		rays[i].origin    = rvs.template unit<3>(2.0);
		rays[i].direction = rvs.template unit<3>();
	}
	// triangle points
	Vec<T,3> p0(0.75,0,0);
	Vec<T,3> p1(0,0.75,0);
	Vec<T,3> p2(0,0,0.75);

	clock_t start = clock();
	for (index_t i = 0; i < iters; i++) {
		trace_tri(p0, p1, p2, rays[i], HIT_FRONT);
	}
	clock_t end = clock();
	
	delete [] rays;

	return (end-start) / (double)CLOCKS_PER_SEC;
}

template <typename T, index_t Bands> double profile_sh(index_t iters) {
    SphericalHarmonics<T,Bands> sh;
    index_t n = std::min(iters, (index_t)1000000);
    Vec<T,3> *vs = new Vec<T,3>[n];
    fill_unit_vec_array<T,3>(vs, n);
    
	clock_t start = clock();
    for (index_t i = 0; i < iters; i++) {
        Vec<T,3> v = vs[i % n];
        sh.eval(v);
    }
	clock_t end = clock();
    
    delete [] vs;
    
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <typename T, index_t Bands> double profile_sh_buffer(index_t iters) {
    SphericalHarmonics<T,Bands> sh;
    SphericalHarmonics<T,Bands> buf;
    index_t n = std::min(iters, (index_t)1000000);
    Vec<T,3> *vs = new Vec<T,3>[n];
    fill_unit_vec_array<T,3>(vs, n);
    T k = 0;
    
	clock_t start = clock();
    for (index_t i = 0; i < iters; i++) {
        Vec<T,3> v = vs[i % n];
        spherical_harmonic_coeff(&buf, v.z, std::atan2(v.y, v.x));
        k = sh.dot(buf);
    }
	clock_t end = clock();
    
    delete [] vs;
    
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <typename T, index_t Bands> 
void test_sh_methods(index_t iters) {
    SphericalHarmonics<T,Bands> sh;
    SphericalHarmonics<T,Bands> buf;
    index_t n = std::min(iters, (index_t)1000000);
    Vec<T,3> *vs = new Vec<T,3>[n];
    fill_unit_vec_array<T,3>(vs, n);
    
    std::fill(sh.coeffs.get(), sh.coeffs.get() + sh.size(), 1);
    
    T diff = 0;
    T a_avg = 0;
    T b_avg = 0;
    
    for (index_t i = 0; i < iters; i++) {
        Vec<T,3> v = vs[i % n];
        spherical_harmonic_coeff(&buf, v.z, std::atan2(v.y, v.x));
        T a = sh.dot(buf);
        T b = sh.eval(v);
        a_avg += a;
        b_avg += b;
        diff  += std::abs(a-b);
    }
    
    std::cout << "avg sh diff: " << diff / iters << std::endl;
    std::cout << "avg buf,std: " << a_avg << "," << b_avg << std::endl;
    
    delete [] vs;
}

template <typename T, index_t N> double profile_path(index_t iters) {
    Path<T,N> p;
    Sampler<T> rvs;
    int n_knots = 50;
    for (index_t i = 0; i < n_knots; i++) {
        p.knots.push_back(Ray<T,N>(rvs.template box<N>(), rvs.template solidball<N>()));
    }
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++) {
        p.eval(getRandom()->rand(0,n_knots));
    }
    clock_t end = clock();
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <typename T, index_t N, index_t Channels, Interpolation Interp> 
double profile_raster(index_t iters) {
	index_t n = 128;
	Vec<index_t,N> dim(n);
    Raster<T,T,N,Channels> image(dim);
    Vec<T,N> *coords = new Vec<T,N>[n];
    Sampler<T> rvs;
    
    for (index_t i = 0; i < n; i++) {
    	coords[i] = rvs.box(Vec<T,N>::zeros, (Vec<T,N>)dim);
    	for (index_t j = 0; j < n; j++) {
    		image.set(Vec<index_t,N>(i,j), rvs.template unit<Channels>());
    	}
    }
    
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++) {
        image.template sample<EDGE_CLAMP,Interp>(coords[i & (n-1)]);
    }
    clock_t end = clock();
    
    delete [] coords;
    
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <typename T> double profile_rotCtr(index_t iters) {
    Sampler<T> rvs;
    AffineTransform<T,3> xf;
    clock_t start = clock();
    for (index_t i = 0; i < iters; i++) {
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
	for (index_t i = 0; i < rows*cols+1; i++) {
		ary[i] = rng->rand(1.0);
	}
	clock_t start = clock();
	for (index_t i = 0; i < iters; i++) {
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
	for (index_t i = 0; i < rows*cols+1; i++) {
		ary[i] = rng->rand(1.0);
	}
	clock_t start = clock();
	for (index_t i = 0; i < iters; i++) {
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
	for (T *p = mtx[0].begin(); p != mtx[0].end(); p++) {
		*p = rng->rand(1.0);
	}
	index_t x = 0;
	clock_t start = clock();
	for (index_t i = 0; i < iters; i++, x = !x) {
		inv(&mtx[!x], mtx[x]);	
	}
	clock_t end = clock();
	
	return (end-start) / (double)CLOCKS_PER_SEC;
}

template <typename T, index_t N> double profile_mtxInverseLU(index_t iters) {
	SimpleMatrix<T,N,N> mtx[2];
	Random *rng = getRandom();
	for (T *p = mtx[0].begin(); p != mtx[0].end(); p++) {
		*p = rng->rand(1.0);
	}
	index_t x = 0;
	clock_t start = clock();
	for (index_t i = 0; i < iters; i++, x = !x) {
		PLUDecomposition<T,N,N> plu = plu_decompose(mtx[x]);
		plu.inverse(&mtx[!x]);	
	}
	clock_t end = clock();
	
	return (end-start) / (double)CLOCKS_PER_SEC;
}

template <typename T, index_t N>
void randomOrient(AffineTransform<T,N> *xf) {
    Sampler<T> smp;
    Vec<T,N> basis[N];
    for (index_t i = 0; i < N; i++) {
        basis[i] = smp.template unit<N>();
    }
    orthonormalize(basis,N);
    SimpleMatrix<T,N,N> rot(basis[0].begin());
    *xf = transformation(rot);
}

template <typename T>
void randomOrient(AffineTransform<T,3> *xf) {
    Sampler<T> smp;
    Quat<T> q = smp.template unit<4>();
    *xf = rotation(q);
}

template <typename T>
void randomOrient(AffineTransform<T,2> *xf) {
    T angle =  M_PI * getRandom()->rand<T>(-1,1);
    *xf = rotation(angle);
}

template <typename T, index_t N>
void randomBox(OrientedRect<T,N> *r) {
    Sampler<T> smp;
    
    randomOrient(&r->xf);
    r->xf *= scale(Vec<T,N>(getRandom()->rand<T>(2,8)));
    r->xf *= translation(smp.template unit<N>());
    Vec<T,N> b0 = smp.box(Vec<T,N>(-1), Vec<T,N>(1));
    Vec<T,N> b1 = smp.box(Vec<T,N>(-1), Vec<T,N>(1));
    r->box = Rect<T,N>::spanningCorners(b0,b1);
}

#ifdef ENABLE_FRUSTUM

template <typename T>
void randomFrustum(Frustum<T,2> *f) {
    Sampler<T> smp;
    f->height = Rect<T,1>::spanningCorners(
                            getRandom()->rand<T>(-5), 
                            getRandom()->rand<T>( 5));
    f->base = Rect<T,1>::spanningCorners(
                            getRandom()->rand<T>(-5), 
                            getRandom()->rand<T>( 5));
    randomOrient(&f->xf);
    f->xf *= scale(Vec<T,2>(getRandom()->rand(2,8)));
    f->xf *= translation(smp.template unit<2>());
}

template <typename T, index_t N>
void randomFrustum(Frustum<T,N> *f) {
    Sampler<T> smp;
    f->height = Rect<T,1>::spanningCorners(
                            getRandom()->rand<T>(-5), 
                            getRandom()->rand<T>( 5));
    Vec<T,N-1> b0 = smp.box(Vec<T,N-1>(-5), Vec<T,N-1>(5));
    Vec<T,N-1> b1 = smp.box(Vec<T,N-1>(-5), Vec<T,N-1>(5));
    f->base = Rect<T,N-1>::spanningCorners(b0, b1);
    randomOrient(&f->xf);
    f->xf *= scale(Vec<T,N>(getRandom()->rand(2,8)));
    f->xf *= translation(smp.template unit<N>());
}

#endif

template <typename T, index_t N> double profile_gjkIntersect(index_t iters) {
    const index_t n = (index_t)std::ceil(std::sqrt(iters));
    const index_t n_corners = 1 << N;
    //Vec<T,N> *rects = new Vec<T,N>[n*n_corners];
    OrientedRect<T,N> *rects = new OrientedRect<T,N>[n];
    bool b = true;
    for (index_t i = 0; i < n; i++) {
        //OrientedRect<T,N> r;
        randomBox(rects+i);
        //r.getCorners(rects + n_corners*i);
    }
    
    Vec<T,N> d;
    index_t i0 = 0;
    index_t i1 = 0;
    clock_t start = clock();
    for (index_t j = 0; j < iters; j++) {
        i0 = (i0 + 1) % n;
        if (i0 == 0) i1 = (i1 + 1) % n;
        //b = b ^ gjk_intersect(rects + n_corners*i0, n_corners, rects + n_corners*i1, n_corners, &d);
        b = b ^ gjk_intersect(rects[i0], rects[i1], &d);
    }
    clock_t end = clock();
    delete [] rects;
    
    return (end-start) / (double)CLOCKS_PER_SEC;
}

template <typename T, index_t N> double profile_OBB_intersect(index_t iters) {
    const index_t n = (index_t)std::ceil(std::sqrt(iters));
    OrientedRect<T,N> *boxes = new OrientedRect<T,N>[n];
    bool b = true;
    for (index_t i = 0; i < n; i++) {
        randomBox(boxes + i);
    }
    
    Vec<T,N> d;
    index_t i0 = 0;
    index_t i1 = 0;
    clock_t start = clock();
    bool SAT = false;
    for (index_t j = 0; j < iters; j++) {
        i0 = (i0 + 1) % n;
        if (i0 == 0) i1 = (i1 + 1) % n;
        SAT = SAT ^ boxes[i0].intersects(boxes[i1]);
    }
    clock_t end = clock();
    
    delete [] boxes;
    
    return (end-start) / (double)CLOCKS_PER_SEC;
}

// this is tautological for N > 3. OBBs use GJK directly!
template <typename T, index_t N> index_t test_gjkIntersect(index_t iters) {
    const index_t n = (index_t)std::ceil(std::sqrt(iters));
    OrientedRect<T,N> *boxes = new OrientedRect<T,N>[n];
    bool b = true;
    for (index_t i = 0; i < n; i++) {
        randomBox(boxes + i);
    }
    const index_t n_corners = 1 << N;
    
    Vec<T,N> d;
    index_t i0 = 0;
    index_t i1 = 0;
    index_t failures = 0;
    index_t positive = 0;
    index_t negative = 0;
    for (index_t j = 0; j < iters; j++) {
        i0 = (i0 + 1) % n;
        if (i0 == 0) i1 = (i1 + 1) % n;
        Vec<T,N> b0[n_corners];
        Vec<T,N> b1[n_corners];
        boxes[i0].getCorners(b0);
        boxes[i1].getCorners(b1);
        bool gjk = gjk_intersect(b0, n_corners, b1, n_corners, &d);
        bool SAT = boxes[i0].intersects(boxes[i1]);
        if (gjk != SAT) failures++;
        if (SAT) positive++;
        else negative++;
    }
    
    delete [] boxes;
    
    std::cout << "gjk " << N << "D failures: " << failures;
    std::cout << " (" << positive << " overlapped " << negative << " disjoint)" << std::endl;
    
    return failures;
}


#ifdef ENABLE_FRUSTUM

template <typename T, index_t N> index_t test_frustumSupport(index_t iters) {
    const index_t n = (index_t)std::ceil(std::sqrt(iters));
    const index_t n_corners = 1 << N;
    
    Sampler<T> smp;
    T err = 0;
    
    index_t failures = 0;
    for (index_t j = 0; j < iters/100; j++) {
        Frustum<T,N> f;
        randomFrustum(&f);
        Vec<T,N> pts[n_corners];
        f.getCorners(pts);
        for (index_t k = 0; k < 100; k++) {
            Vec<T,N> d = smp.template unit<N>();
            Vec<T,N> p0 = f.convexSupport(d);
            Vec<T,N> p1 = pts[0];
            T dot0 = p1.dot(d);
            // brute force find support point
            for (index_t c = 1; c < n_corners; c++) {
                T dot1 = pts[c].dot(d);
                if (dot1 > dot0) {
                    dot0 = dot1;
                    p1 = pts[c];
                }
            }
            if (p0 != p1) {
                err += p1.dist(p0);
                failures++;
            }
        }
    }
    
    std::cout << "frustum support " << N << "D failures: " << failures;
    std::cout << " (" << (100 * failures / (double)iters) << "%)";
    std::cout << " average mismatch: " << (err / (double)iters) << std::endl;
    
    return failures;
}

#endif

template <typename T, index_t N> void test_mtxInverse(index_t iters) {
	SimpleMatrix<T,N,N> mtx[2];
	Random *rng = getRandom();
	for (T *p = mtx[0].begin(); p != mtx[0].end(); p++) {
		*p = rng->rand(1.0);
	}
	std::cout << mtx[0];
	index_t x = 0;
	for (index_t i = 0; i < iters; i++, x = !x) {
		inv(&mtx[!x], mtx[x]);
		std::cout << mtx[!x];
	}
}

template <typename M> void test_mtxRegionCopy(index_t rows, index_t cols) {
	M m(rows,cols);
	typedef typename M::elem_t T;
	const index_t subn = (rows-2)*(cols-2);
	T *v = new T[subn];
	
	for (index_t i = 0; i < subn; i++) {
		v[i] = 11*i;
	}
	
	std::fill(m.begin(), m.end(), 1);
	std::copy(v, v+subn, m.region_begin(1,1,rows-1, cols-1));
	cout << m;
}

template <typename T, index_t N> void test_OOBB() {
    OrientedRect<T,N> r1;
    OrientedRect<T,N> r2;
    cout << "box1 intersects box2: " << r1.intersects(r2) << endl;
}

template <typename T, index_t M, index_t N> void test_mtxArithmetic(index_t rows, index_t cols) {
	SimpleMatrix<T,M,N> a(rows,cols);
    SimpleMatrix<T,M,N> b(rows,cols);
    SimpleMatrix<T,M,N == 0 ? 0 : N+1> bogus(rows,cols+1);
    
    cout << "addTest:" << endl;
    cout << a + b << endl;
    cout << a * 3 << endl;
    cout << 3 * a << endl;
    //correctly does not compile:
    //cout << a + bogus << endl;
}

template <index_t N> void test_permuteMatrix() {
	PermutationMatrix<N> p;
    PermutationMatrix<N> p1;
	Vec<double,N> v;
	Random *rng = getRandom();
	for (index_t i = 0; i < N; i++) {
		v[i] = i;
	}
	for (index_t i = 0; i < N*2; i++) {
		p.swap_rows(rng->rand(N), rng->rand(N));
        p1.swap_rows(rng->rand(N), rng->rand(N));
	}
	SimpleMatrix<double,N,N> m;
	SimpleMatrix<double,N,N> minv;
    SimpleMatrix<double,N,N> m_p1;
	std::copy(p.begin(), p.end(),  m.begin());
    std::copy(p1.begin(), p1.end(),  m_p1.begin());
	inv(&minv, p);
	cout << "v: " << v << endl;
	cout << "P: " << endl << p << endl;
	cout << "P*v: " << (p * v) << endl;
	cout << "M: " << endl << m << endl;
	cout << "M*v: " << (m * v) << endl;
    cout << "P * P1: " << (p * p1) << endl;
    cout << "P * P1 correct? " << ((p * p1) == (m * m_p1)) << endl;
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

void profile(std::string name, double (*fnc)(index_t), index_t iterations) {
    double t = fnc(iterations);
    std::cout << name << " (" << iterations << " iters): " << t << " secs; " << iterations/t << " ops/sec"  << std::endl;
}

void test_matrix_asplode() {
    SimpleMatrix<double,4,4> mtx;
    SimpleMatrix<double,2,5> other;
    // correctly, won't compile:
    // mtx * other;
}

void test_dual() {
    typedef Dual<double> duald;
    cout << "sin(pi): " << sin(duald(PI,1)) << endl;
    cout << "cos(pi): " << cos(duald(PI,1)) << endl;
    cout << "tan(pi): " << tan(duald(PI,1)) << endl;
    cout << "2x + 1, @x=3: " << (duald(3,1) * 2 + 1) << endl;
    cout << "exp(1): " << exp(duald(1,1)) << endl;
    cout << "pow(x,x) @x=2: " << pow(duald(2,1), duald(2,1)) << endl;
    cout << "pow(x=2,2): " << pow(duald(2,1), 2) << endl;
    cout << endl;
}

using geom::operator<<;


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
    test_mtxArithmetic<double,4,3>(4,3);
    test_dual();
    test_OOBB<double,3>();
    test_gjkIntersect<double,2>(1000000);
    test_gjkIntersect<double,3>(1000000);
    
#ifdef ENABLE_FRUSTUM
    test_frustumSupport<double,2>(1000000);
    test_frustumSupport<double,3>(1000000);
    test_frustumSupport<double,4>(1000000);
    std::cout << std::endl;
#endif
    
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
    
    //profile("1d perlin dual", profile_perlin_grad<double,1>, iters/100);
    profile("2d perlin dual", profile_perlin_grad<double,2>, iters/100);
    profile("3d perlin dual", profile_perlin_grad<double,3>, iters/100);
    profile("4d perlin dual", profile_perlin_grad<double,4>, iters/100);
    std::cout << std::endl;
    
    //profile("1d perlin grad", profile_perlin_hard_grad<double,1>, iters/100);
    profile("2d perlin grad", profile_perlin_hard_grad<double,2>, iters/100);
    profile("3d perlin grad", profile_perlin_hard_grad<double,3>, iters/100);
    profile("4d perlin grad", profile_perlin_hard_grad<double,4>, iters/100);
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
    
    profile("sh  3f band", profile_sh<float, 3>,  iters/10);
    profile("sh  8f band", profile_sh<float, 8>,  iters/50);
    profile("sh 16f band", profile_sh<float, 16>, iters/100);
    profile("sh  3d band", profile_sh<double, 3>,  iters/10);
    profile("sh  8d band", profile_sh<double, 8>,  iters/50);
    profile("sh 16d band", profile_sh<double, 16>, iters/100);
    std::cout << std::endl;
    
    profile("2f->3f linear img sample", profile_raster<float, 2, 3, INTERP_LINEAR>, iters/100);
    profile("3f->3f linear img sample", profile_raster<float, 3, 3, INTERP_LINEAR>, iters/100);
    profile("2f->3f cubic img sample",  profile_raster<float, 2, 3, INTERP_CUBIC>,  iters/100);
    profile("3f->3f cubic img sample",  profile_raster<float, 3, 3, INTERP_CUBIC>,  iters/100);
    std::cout << std::endl;
    
    profile("5x5 mtxf region copy",    profile_mtxRegionCopy<SimpleMatrix<float,0,0> >,    iters);
    profile("5x5 mtxd region copy",    profile_mtxRegionCopy<SimpleMatrix<double,0,0> >,   iters);
    profile("5x5 sparsef region copy", profile_mtxRegionCopy<SparseMatrix<float> >,  iters/100);
    profile("5x5 sparsed region copy", profile_mtxRegionCopy<SparseMatrix<double> >, iters/100);
    profile("5x5 diagf region copy",   profile_mtxRegionCopy<DiagMatrix<float,5,5> >, iters);
    std::cout << std::endl;

    profile("5x5 mtxf copy",    profile_mtxCopy<SimpleMatrix<float,0,0> >, iters/100);
    profile("5x5 mtxd copy",    profile_mtxCopy<SimpleMatrix<double,0,0> >, iters/100);
    profile("5x5 sparsef copy", profile_mtxCopy<SparseMatrix<float> >, iters/100);
    profile("5x5 sparsed copy", profile_mtxCopy<SparseMatrix<double> >, iters/100);
    profile("5x5 diagf copy",   profile_mtxCopy<DiagMatrix<float,5,5> >, iters/100);
    std::cout << std::endl;
    
    profile("2x2 mtxf inv", profile_mtxInverse<float,2>,  iters/10);
    profile("3x3 mtxf inv", profile_mtxInverse<float,3>,  iters/10);
    profile("4x4 mtxf inv", profile_mtxInverse<float,4>,  iters/10);
    profile("2x2 mtxd inv", profile_mtxInverse<double,2>, iters/10);
    profile("3x3 mtxd inv", profile_mtxInverse<double,3>, iters/10);
    profile("4x4 mtxd inv", profile_mtxInverse<double,4>, iters/10);
    profile("8x8 mtxd inv", profile_mtxInverse<double,8>, iters/10);
    std::cout << std::endl;
    
    profile("rotationd matrix", profile_rotCtr<double>, iters/10);
    profile("rotationf matrix", profile_rotCtr<float>, iters/10);
    std::cout << std::endl;
    
    profile("gjk float 2d",  profile_gjkIntersect<float,2>,  1000000);
    profile("gjk double 2d", profile_gjkIntersect<double,2>, 1000000);
    profile("gjk float 3d",  profile_gjkIntersect<float,3>,  1000000);
    profile("gjk double 3d", profile_gjkIntersect<double,3>, 1000000);
    profile("gjk float 4d",  profile_gjkIntersect<float,4>,  100000);
    profile("gjk double 4d", profile_gjkIntersect<double,4>, 100000);
    profile("gjk double 5d", profile_gjkIntersect<double,5>, 100000);
    std::cout << std::endl;
    
    profile("SAT 2D OBB double", profile_OBB_intersect<double,2>, 1000000);
    profile("SAT 3D OBB double", profile_OBB_intersect<double,3>, 1000000);
    std::cout << std::endl;

    std::cout << "sizeof vec2f: " << sizeof(Vec<float,2>) << std::endl;
    std::cout << "sizeof vec3f: " << sizeof(Vec<float,3>) << std::endl;
    std::cout << "sizeof vec4f: " << sizeof(Vec<float,4>) << std::endl;
    std::cout << "sizeof vec2d: " << sizeof(Vec<double,2>) << std::endl;
    std::cout << "sizeof vec3d: " << sizeof(Vec<double,3>) << std::endl;
    std::cout << "sizeof vec4d: " << sizeof(Vec<double,4>) << std::endl;
    
    try {
        test_matrix_asplode();
    } catch (GeomException &ex) {
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
