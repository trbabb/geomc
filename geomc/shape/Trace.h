/*
 * Trace.h
 *
 *  Created on: Mar 6, 2013
 *      Author: tbabb
 */

#ifndef TRACE_H_
#define TRACE_H_

#include <geomc/shape/shapedetail/Hit.h>
#include <geomc/shape/shapedetail/TraceImpl.h>
#include <geomc/shape/ShapeTypes.h>
#include <geomc/shape/Plane.h>
#include <geomc/linalg/LUDecomp.h>
 
// TODO: planar figures, to plug in with trace_planar_quad

namespace geom {


/** @addtogroup shape 
 *  @{ 
 */

/********************************
 * Static tracing functions     *
 ********************************/
 

/**
 * Compute the constructive solid geometry intersection of two convex shapes 
 * from their near and far hits.
 *
 * The given near and far hits may be distances or ray parameters, and
 * may be negative.
 *
 * @param a_lo Near hit of object A.
 * @param a_hi Far hit of object A.
 * @param b_lo Near hit of object B.
 * @param b_hi Far hit of object B.
 * @param s Ray parameter. (return variable)
 * @param hit_a Will be set to `true` if object A's surface generated the 
 * hit, otherwise `false`. (return variable)
 * @param hit_hear Will be set to `true` if the near surface (lower ray parameter)
 * generated the hit, otherwise `false`. (return variable)
 * @param side Return variable for the surface side (HIT_FRONT or HIT_BACK) that was
 * hit; must be initialized with the permitted hit side flags.
 */
template <typename T, index_t N>
inline bool ray_csg_intersect(
        T a_lo, T a_hi,
        T b_lo, T b_hi,
        T* s,
        bool* hit_a,
        bool* hit_near,
        HitSide* side) {
    
    // csg intersect: take the farthest near and nearest far
    bool a_near = a_lo > b_lo;
    bool a_far  = a_hi < b_hi;
    T s_lo      = a_near ? a_lo : b_lo;
    T s_hi      = a_far  ? a_hi : b_hi;
    
    // choose near or far hit according to hit side
    if (s_lo > 0 && (*side & HIT_FRONT)) {
        *hit_near = true;
        *side     = HIT_FRONT;
        *hit_a    = a_near;
    } else if (s_hi > 0 && (*side & HIT_BACK)) {
        *hit_near = false;
        *side     = HIT_BACK;
        *hit_a    = a_far;
    } else {
        // side test failed; no hit
        return false;
    }
    return true;
}

/**
 * Compute the constructive solid geometry subtraction (A-B) of two 
 * convex shapes from their near and far hits. 
 *
 * The given near and far hits may be distances or ray parameters, and
 * may be negative.
 *
 * @param a_lo Near hit of object A.
 * @param a_hi Far hit of object A.
 * @param b_lo Near hit of subtracted object B.
 * @param b_hi Far hit of subtracted object B.
 * @param s Ray parameter. (return variable)
 * @param hit_a Will be set to `true` if object A's surface generated the 
 * hit, otherwise `false`. (return variable)
 * @param hit_hear Will be set to `true` if the near surface of the constructed
 * object produced the hit; `false` otherwise. (return variable)
 * @param side Return variable for the surface side (HIT_FRONT or HIT_BACK) that was
 * hit; must be initialized with the permitted hit side flags.
 */
template <typename T, index_t N>
inline bool ray_csg_subtract(
        T a_lo, T a_hi,
        T b_lo, T b_hi,
        T* s,
        bool* hit_a,
        bool* hit_near,
        HitSide* side) {
    // xxx: I don't think this is right.
    // consider the concave part of a subtracted sphere.
    
    // csg subtract
    bool a_near = a_lo < b_hi;
    bool a_far  = a_hi < b_lo;
    T s_lo      = a_near ? a_lo : b_hi;
    T s_hi      = a_far  ? a_hi : b_lo;
    
    // choose near or far hit according to hit side
    if (s_lo > s_hi) {
        // empty interval; no hit
        return false;
    } else if (s_lo > 0 && (*side & HIT_FRONT)) {
        *hit_near = true;
        *side     = HIT_FRONT;
        *hit_a    = a_near;
    } else if (s_hi > 0 && (*side & HIT_BACK)) {
        *hit_near = false;
        *side     = HIT_BACK;
        *hit_a    = a_far;
    } else {
        // side test failed; no hit
        return false;
    }
    return true;
}


/**
 * Ray trace a simplex (e.g. triangle, tetrahedron, etc.), returning
 * the surface parameters of the hit point, given in coordinates of the basis 
 * spanned by the edges radiating from `verts[0]`.
 * 
 * @param verts An array of the simplex's N vertices.
 * @param ray The ray to intersect with the simplex.
 * @param uv Buffer for the surface coordinates of the hitpoint. (return variable)
 * @param s The ray parameter of the hit point. (return variable)
 * @return `true` if the ray hit the simplex.
 */
template <typename T, index_t N>
bool trace_simplex(const Vec<T,N> verts[N], const Ray<T,N>& ray, Vec<T,N-1>* uv, T* s) {
    /* Solve for u,v coordinates along the edges radiating from `verts[0]`:
    A + u(B-A) + v(C-A)      = o + sV
        u(B-A) + v(C-A) - sV = o - A
    Mx = o - A
    */
    
    // todo: cramer's rule might be used in 3d for faster linear solve. 
    //       (note: is unstable).
    // todo: can formulate row-major for speed?
    
    // populate matrix
    SimpleMatrix<T,N,N> m;
    for (index_t col = 1; col < N; ++col) {
        Vec<T,N> v = verts[col] - verts[0];
        for (index_t row = 0; row < N; ++row) {
            m.set(row, col, v[row]);
        }
    }
    for (index_t row = 0; row < N; ++row) {
        m.set(row, N-1, -ray.direction[row]);
    }
    
    // linear solve
    Vec<T,N> x;
    Vec<T,N> b = ray.origin - verts[0];
    if (N < 5) {
        SimpleMatrix<T,N,N> m_inv;
        if (!inv(&m_inv, m)) return false;
        x = m_inv * b;
    } else {
        if (!linearSolve(m.begin(), N, x.begin(), b.begin())) return false;
    }
    
    // inside simplex?
    T sum = 0;
    for (index_t i = 0; i < N - 1; ++i) {
        if (x[i] < 0) return false;
        sum += x[i];
    }
    if (sum > 1) return false;
    
    // output result
    *s  = x[N];
    *uv = x.template resized<N-1>();
    return true;
}


template <typename T, index_t N>
inline Hit<T,3> trace_tri(Vec<T,3> p[3], const Ray<T,3> &r, HitSide sides=HIT_FRONT){
    return trace_tri(p[0], p[1], p[2], r, sides);
}


// todo: this is undoubtedly slower than it has to be
template <typename T>
Hit<T,3> trace_tri(const Vec<T,3> &p0, 
                   const Vec<T,3> &p1, 
                   const Vec<T,3> &p2, 
                   const Ray<T,3> &r, 
                   HitSide sides=HIT_FRONT){
    Hit<T,3> hit(r,sides);
    
    Vec<T,3> dP[2]   = {p1 - p0, p2 - p0};
    Plane<T,3> plane = Plane<T,3>::from_basis(dP, p0);
    T s;
    
    if (!detail::_ImplTracePlane(&s, plane, r, &sides)) {
        return hit;
    }
    
    Vec<T,3> p = r.atMultiple(s);
    
    // dPdu, dPdv, and N form a 3d basis.
    // thus we can solve:
    //    u*dPdu + v*dPdv + t*N = P
    // t should essentially be zero; we can ignore it.
    //     xxx: ^ that reeks of stupid. can we solve this in one go?
    
    // 3d matrix inversion is faster than you think. 
    // also, gcc tends to unroll those copy ops.
    SimpleMatrix<T,3,3> mtx;
    std::copy(dP[0].begin(), dP[0].end(), mtx.col(0));
    std::copy(dP[1].begin(), dP[1].end(), mtx.col(1));
    std::copy(plane.normal.begin(), plane.normal.end(), mtx.col(2));
    inv(&mtx, mtx);
    Vec<T,3> uvn = mtx * (p-p0);
    
    T u = uvn.x;
    T v = uvn.y;
    
    if (u < 0 || v < 0 || u + v > 1) {
         //intersection outside of triangle
        return hit;
    }
    
    hit.p = p;
    hit.n = plane.normal;
    hit.s = s;
    hit.side = sides; //set by ImplTracePlane
    hit.hit = true;
    
    return hit;
}


template <typename T>
Hit<T,3> trace_planar_quad(
        const Vec<T,3> &dPdu, 
        const Vec<T,3> &dPdv, 
        const Vec<T,3> &origin, 
        const Ray<T,3> &ray, 
        HitSide sides) {
    Hit<T,3> hit(ray,sides);
    Vec<T,3> dP[2] = {dPdu, dPdv};
    Plane<T,3> pl = Plane<T,3>::from_basis(dP, origin);
    T s;
    
    if (!detail::_ImplTracePlane(&s, pl, ray, &sides)){
        // miss
        return hit;
    }
    
    Vec<T,3> p = ray.atMultiple(s);
    
    // solve for u and v
    SimpleMatrix<T,3,3> mtx;
    std::copy(dPdu.begin(), dPdu.end(), mtx.col(0));
    std::copy(dPdv.begin(), dPdv.end(), mtx.col(1));
    std::copy(pl.normal.begin(), pl.normal.end(), mtx.col(2));
    inv(&mtx, mtx);
    Vec<T,3> uvn = mtx * (p-origin);
    
    T u = uvn.x;
    T v = uvn.y;
    
    if (u < 0 or u > 1 or v < 0 or v > 1){
        // missed quad area
        return hit;
    }
    
    hit.p = p;
    hit.n = pl.normal;
    hit.s = s;
    hit.side = sides; // set by ImplTracePlane
    hit.hit = true;
    
    return hit;
}

    
template <typename T, index_t N>
inline T trace_plane(const Vec<T,N> plane_pt, const Vec<T,N> plane_n, const Ray<T,N> &r) {
    Vec<T,N> w = r.origin - plane_pt;
    T z = w.dot(plane_n);
    T s = -z / r.direction.dot(plane_n);
    return s;
}


template <typename T, index_t N>
inline T trace_plane(const Vec<T,N> plane_pt, const Vec<T,N> plane_n, const Ray<T,N> &r, T* dot) {
    Vec<T,N> w = r.origin - plane_pt;
    T z  = w.dot(plane_n);
    T s  = -z / r.direction.dot(plane_n);
    *dot = z; // which side of the plane? + --> side of normal
    return s;
}

template <typename T, index_t N>
bool trace_sphere(const Vec<T,N>& center, T radius, const Ray<T,N> &ray, T* s) {
    T r2 = radius * radius;
    // if inside, we are guaranteed to hit the back. return miss if not tracing backface.
    const Vec<T,N>& dir = ray.direction; //for shorthand
    Vec<T,N> x0 = center - ray.origin;
    // solve for s such that ||s * ray - ctr|| == radius 
    T a = dir.dot(dir);
    T b = -2 * dir.dot(x0);
    T c = x0.dot(x0) - r2;
    T roots[2];
    if (quadratic_solve(roots, a, b, c)) {
        if (roots[1] < roots[0]) {
            std::swap(roots[0], roots[1]);
        }
        T positive_root = (roots[0] >= 0) ? roots[0] : roots[1];
        if (positive_root >= 0) {
            *s = positive_root;
            return true;
        }
    }
    // no intersection; return miss.
    return false;
}


/// @}   // addtogroup shape


} // namespace geom

#endif /* TRACE_H_ */
