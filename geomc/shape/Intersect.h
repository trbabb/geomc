/* 
 * File:   Intersect.h
 * Author: tbabb
 *
 * Created on October 14, 2014, 11:30 PM
 */

#ifndef INTERSECT_H
#define	INTERSECT_H

#include <limits.h>
#include <geomc/linalg/Vec.h>
#include <geomc/linalg/Orthogonal.h>

#include "Plane.h"

namespace geom {
    
namespace detail {
    
    template <typename T, index_t N>
    struct Simplex {
        Simplex():n(0) {}
        
        Vec<T,N> pts[N+1];
        index_t  n;
        
        inline void insert(const Vec<T,N> &p) {
            if (n <= N) pts[n++] = p;
        }
        
        inline void remove(index_t i) {
            for (index_t c = i; c < n - 1; c++) {
                pts[c] = pts[c+1];
            }
            n -= 1;
        }
    };
    
    // todo: can we cache partial evaluations of LUP()?
    // todo: can we pass in/save the null space?
    
    // todo: profile this. compare to separating axis.
    //       verify whether final projection space swap actually helps or not.
    
    // return true if `s` contains origin
    // otherwise, put the direction toward the origin in `d`.
    template <typename T, index_t N>
    bool gjk_simplex_near_origin(Simplex<T,N> *s, Vec<T,N> *d) {
        // let the most recently added vertex in `s` be called `A`. 
        // we need not check simplexes and sub-simplexes that don't involve
        // A, because they were checked in previous iterations of GJK, before 
        // A was added. 
        Vec<T,N> A = s->pts[s->n - 1];
        
        // find the null space of our simplex (the normals to sub-simplexes will
        // be orthogonal to it). the high part of this array holds the basis of 
        // the simplex we're testing:
        Vec<T,N>  simplex_null_basis[N]; 
        Vec<T,N> *simplex_basis = simplex_null_basis + (N - (s->n - 1));
        if (s->n > 1 and s->n < N + 1) {
            for (int i = 0; i < s->n - 1; i++) {
                simplex_basis[i] = s->pts[i] - A;
            }
            nullspace(simplex_basis, s->n - 1, simplex_null_basis);
        }
        
        // find the boundary simplex.
        bool removed_vertex = false;
        for (index_t i = 0; i < s->n - 1; ) {
            // pick a candidate vertex to exclude from the simplex, creating a 
            // simplex with one fewer dimension. is this sub-simplex a possible
            // boundary simplex? if so, we discard the vertex and recurse.
            Vec<T,N> excluded = s->pts[i];
            Vec<T,N> normal;
            
            if (s->n > 2) {
                // compute the normal, which is orthogonal to the basis of 
                // s_candidate and also to the null space of s_parent. recall 
                // that simplex_null_basis holds both!
                index_t ctr = 0;
                for (index_t j = 0; j < s->n - 1; j++) {
                    if (j == i) continue;
                    simplex_basis[ctr] = s->pts[j] - A;
                    ctr += 1;
                }
                normal = orthogonal(simplex_null_basis);
            } else {
                // our simplex is a line, and we are testing the endpoint.
                // (the above would also do the right thing, but more slowly).
                normal = A - excluded;
            }
            
            // does the normal we picked point "inside"? if so, flip it
            if (normal.dot(excluded - A) > 0) normal = -normal;
            
            // is the origin in the direction of the normal we picked?
            if (A.dot(normal) < 0) {
                // the current sub-simplex (i.e. the one excluding the candidate 
                // vertex) must contain or itself be the boundary simplex; therefore
                // we can safely remove the candidate vertex.
                s->remove(i);
                // the dimensionality of our boundary simplex is reduced;
                // we have a new axis for our null space:
                simplex_basis[0] = normal;
                simplex_basis++;
                removed_vertex = true;
            } else {
                // we're keeping the candidate vertex. Move on to the next.
                i += 1;
            }
        }
        
        index_t null_dim = N - (s->n - 1);
        
        // find the direction to the origin
        // i.e. project the direction to the origin onto the null space of our simplex.
        // the final `else` would be sufficient by itself, but in some cases
        // we can be clever and save some redundant calculations.
        if (null_dim == N) {
            // single point case; projection is trivial.
            *d = -A;
        } else if (null_dim == 1) {
            // the simplex spans a hyperplane; we already have the normal.
            *d = simplex_null_basis[0];
            // if we created the plane by removing a vertex, we already know the 
            // normal is facing toward the origin. otherwise, we need to check 
            // and flip it if necessary.
            if (!removed_vertex and d->dot(A) > 0) *d *= -1;
        } else if (null_dim == 0) {
            // simplex forms a complete basis; no vertexes were eliminated.
            // in other words, no plane test failed. The origin is inside
            // this simplex!
            return true;
        } else {
            // explicitly project the direction to the origin onto the null space
            // of the boundary simplex.
            *d = Vec<T,N>::zeros;
            
            // either project onto the null basis directly, or subtract a vector
            // orthogonal to the null basis from -A, according to whichever 
            // projection is simpler.
            bool swap_proj = false;
            index_t proj_n = null_dim;
            Vec<T,N> *proj_basis = simplex_null_basis;
            
            if ( null_dim > s->n - 1 ) {
                proj_basis = simplex_basis;
                proj_n = s->n - 1;
                swap_proj = true;
                for (index_t i = 0; i < s->n - 1; i++) {
                    proj_basis[i] = s->pts[i] - A;
                }
            }
            
            // gram-schmidt orthonormalization.
            for (index_t b = 0; b < proj_n; b++) {
                for (index_t i = 0; i < b; i++) {
                    proj_basis[b] -= proj_basis[b].projectOn(proj_basis[i]);
                }
                *d += -A.projectOn(proj_basis[b]);
            }
            
            if (swap_proj) *d = -A - *d;
        }
        return false;
    }
    
    
    // find the point having the largest dot product with direction `d`.
    template <typename T, index_t N>
    Vec<T,N> gjk_support(const Vec<T,N> *pts, index_t n, Vec<T,N> d) {
        T largest_dot = std::numeric_limits<T>::lowest();
        Vec<T,N> pt;
        for (index_t i = 0; i < n; i++) {
            T a = d.dot(pts[i]);
            if (a > largest_dot) {
                largest_dot = a;
                pt = pts[i];
            }
        }
        return pt;
    }

} // namespace detail


// todo: template GJK entirely over the input shape type.

// Use the Gilbert-Johnson-Keerthi algorithm to test for overlap between
// two convex shapes, and return the axis of overlap.
template <typename T, index_t N>
bool gjk_intersect(const Vec<T,N> *pts_a, index_t n_a, 
                   const Vec<T,N> *pts_b, index_t n_b, 
                   Vec<T,N> *overlap_axis) {
    // choose an arbitrary initial direction.
    Vec<T,N> initial;
    if (overlap_axis and *overlap_axis != Vec<T,N>::zeros) 
        initial = *overlap_axis;
    else 
        initial[0] = 1;
    
    Vec<T,N> a = detail::gjk_support<T,N>(pts_a, n_a,  initial) - 
                 detail::gjk_support<T,N>(pts_b, n_b, -initial);
    Vec<T,N> d = -a;
    detail::Simplex<T,N> s;
    s.insert(a);
    
    while (true) {
        a = detail::gjk_support<T,N>(pts_a, n_a,  d) - 
            detail::gjk_support<T,N>(pts_b, n_b, -d);
        if (a.dot(d) < 0) return false;
        s.insert(a);
        if (detail::gjk_simplex_near_origin(&s, &d)) {
            if (overlap_axis) *overlap_axis = d;
            return true;
        }
    }
}

} // end namespace geom

#endif	/* INTERSECT_H */

