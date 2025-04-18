#pragma once

#include <geomc/shape/Shape.h>
#include <geomc/linalg/Vec.h>
#include <geomc/linalg/Similarity.h>
#include <geomc/linalg/Orthogonal.h>
#include <geomc/linalg/Matrix.h>
#include <geomc/shape/shapedetail/SimplexProject.h>

// xxx: contains() is likely wrong
//   - we do need to solve the last coordinate; there are as many coordiantes as verts
// todo: wbn if barycentric() / unmap() projected to simplex subspace and gave coords
//    anyway; then contains() could call this
// todo: boundary measure

// todo: orthogonal projection to a subspace is a very general linalg operation,
//       and projection_contains() makes use of it. The operation that does this is
//       the pseudoinverse, which can be computed efficiently using a QR decomposition.
//       Probably the best way to do QR is with Householder transformations; the algorithm
//       can be found complete here: http://people.inf.ethz.ch/gander/papers/qrneu.pdf
//       
//       the above should probably be exposed as orthogonal_project(basis[], n, P) in
//       Orthogonal.h, and invoked here by projection_contains(). Also consider adding
//       subspace_projection(bases[], n) -> Mtx, which can amortize the QR bidness.
//
//       right now we are running two LUP decompositions, and this is likely
//       inefficient/unstable.
//
//       nb: this might not be useful in lower dimensions

// idea:
//       for the full simplex, solving barycentric coords will give you an implicit
//       face-orthogonal basis: the weight of the vertex "missing" from a face must be 1
//       at the vertex, and 0 across the entire opposing face— thus that coordinate must
//       define isosurfaces parallel to the opposing face!
//       can we solve once and then inspect the solved coordinate values to figure out
//       sub-simplex containment?
//       
//       the trick is to do this for incomplete simplexes. what constraints can we impose
//       to make the linear system square again? (guess: certain dot products = 0)
//       > might also be possible using an O(N) process analogous to gram-schmidt
//       
//       also of concern: does this work for arbitrary sub-simplexes? i.e. do you
//       orthogonally project to an edge by coordinate clipping, or do you have to
//       re-solve?
//       > it might be. the trick is to express scaled barycentric coordinates in
//         terms of un-scaled coordinates.

// todo: another possible tool is to project the simplex to a lower dimension
//       (by dropping coordinates) and barycentric-solving in that lower space.
//       this amounts to projecting the simplex to an axial (hyper-)plane or an axis.
//       the weightings of the verticies will not be different if you add back in
//       the missing coordinates, so the barycentric solution from this lower space
//       is still valid. (This also explains how keep the linear system square).
//       generally you want to project to the axial plane wherein the simplex
//       has the largest projected area. overall this might be more numerically robust.
//       note that you have to find the coordinates of the projected point `p`
//       in the axial plane; this is sort of a double-projection.
//  edit: actually, this doesn't make sense. this doesn't do an orthogonal projection

// todo: handle degeneracy. GJK can make use of it.
// todo: performance check old GJK vs new
//   perf issues:
//     - faster when we know it's always the origin we're testing
//       (can we transform the problem into one where we are checking the origin?)
//     - faster if we know to avoid checking the final face
//       (i.e. the one opposing the excluded vtx)
//     - faster if we indirect the verts
// todo: performance check projection_contains() vs. project() == this


namespace geom {

// fwd decl
template <typename T, index_t N>
bool trace_simplex(const Vec<T,N> verts[N], const Ray<T,N>& ray, Vec<T,N-1>* uv, T* s);


/**
 * @ingroup shape
 * @brief A simplex in N dimensions (e.g. a tetrahedron, triangle, line, point).
 *
 * The simplex may contain up to `N + 1` points, in which case it encloses a 
 * volume, if it is not degenerate. If it contains fewer points, then it
 * spans a subspace of the space in which it is embedded.
 */
template <typename T, index_t N>
class Simplex: public Dimensional<T,N> {
public:

    /// Vertices of this simplex.
    Vec<T,N> pts[N + 1];
    /// Number of vertices in this simplex.
    index_t n;
    
    
    /// Construct an empty simplex, with no vertices.
    constexpr Simplex():n(0) {}
    
    /// Construct an `n`-cornered simplex with vertices at `verts`.
    Simplex(const Vec<T,N>* verts, index_t n):n(n) {
        std::copy(verts, verts + n, pts);
    }
    
    Simplex(std::initializer_list<Vec<T,N>> verts):
            n(std::min<size_t>(verts.size(), N + 1))
    {
        std::copy(verts.begin(), verts.begin() + n, pts);
    }
    
    /// Construct a simplex with its vertices on the unit axes.
    static constexpr Simplex standard_simplex() {
        Simplex s;
        for (index_t i = 1; i < N + 1; ++i) {
            s.pts[i][i - 1] = 1;
        }
        s.n = N + 1;
        return s;
    }
    
    /** 
     * @brief Construct a simplex with all edges of unit length and
     * one edge along the x-axis.
     * 
     * `n` is the number of vertices in the simplex.
     */
    static constexpr Simplex regular_simplex(index_t n=N+1) {
        Simplex s;
        Vec<T,N> c = s.pts[0] = {};
        n = std::min(n, N + 1);
        for (index_t i = 1; i < n; ++i) {
            // compute the barycenter of the "base" of the simplex
            Vec<T,N> c_i = c / i;
            // remaining coordinate is raised above the center of the base by `h`
            T h = std::sqrt(1 - c_i.mag2());
            s.pts[i] = {c_i, h};
            c += s.pts[i];
        }
        s.n = n;
        return s;
    }
    
    static constexpr bool admits_cusps() { return true; }
    
    /// Get the `i`th vertex in this simplex.
    Vec<T,N>& operator[](index_t i) {
        return pts[i];
    }
    
    /// Get the `i`th vertex in this simplex.
    Vec<T,N> operator[](index_t i) const {
        return pts[i];
    }
    
    /**
     * @brief Simplex equality check.
     * 
     * Simplexes are equal if they have identical vertices in identical order.
     */
    bool operator==(const Simplex<T,N>& other) const {
        // overridden because there may be garbage in the 
        // unused vertices which shouldn't count as a difference.
        if (other.n != n) return false;
        for (index_t i = 0; i < n; ++i) {
            if (other.pts[i] != pts[i]) return false;
        }
        return true;
    }
    
    /// Simplex inequality check.
    inline bool operator!=(const Simplex<T,N>& other) const {
        return not ((*this) == other);
    }
    
    /**
     * @brief Return a copy of this simplex with an additional vertex at `p`.
     *
     * If this simplex already has `N + 1` points, a copy of this simplex is returned.
     */
    Simplex<T,N> operator|(const Vec<T,N>& p) const {
        Simplex<T,N> s = *this;
        s.insert(p);
        return s;
    }

    /**
     * @brief Extend this simplex to include `p` by adding `p` as a vertex.
     *
     * Alias for `insert(p)`.
     * 
     * If this Simplex already has `N+1` vertices, then this operator has no effect.
     */
    Simplex<T,N>& operator|=(const Vec<T,N>& p) {
        insert(p);
        return *this;
    }
    
    /**
     * @brief Apply a transformation to the points of this Simplex.
     */
    Simplex<T,N>& operator*=(const AffineTransform<T,N>& xf) {
        for (index_t i = 0; i < n; ++i) {
            pts[i] = xf * pts[i];
        }
        return *this;
    }
    
    /**
     * @brief Apply an inverse transformation to the points of this Simplex.
     */
    Simplex<T,N>& operator/=(const AffineTransform<T,N>& xf) {
        for (index_t i = 0; i < n; ++i) {
            pts[i] = pts[i] / xf;
        }
        return *this;
    }
    
    /**
     * @brief Extend this simplex to include `p` by adding `p` as a vertex.
     * 
     * If this Simplex already has `N+1` vertices, then this function has no effect.
     */
    void insert(Vec<T,N> p) {
        if (n < N + 1) {
            pts[n] = p;
            n += 1;
        }
    }
    
    /// Compute the barycenter of this simplex.
    Vec<T,N> barycenter() const {
        Vec<T,N> c = {};
        for (index_t i = 0; i < n; ++i) {
            c += pts[i];
        }
        return c / n;
    }
    
    /**
     * @brief Ray-shape intersection.
     * 
     * Return the interval over which the ray intersects the simplex.
     * 
     * If the simplex is N-1 dimensional (i.e., planar, like a triangle in 3D), then an
     * intersecting interval will have both of its limits equal to the single value of `s`
     * for which the ray `o + s*v` intersects the simplex.
     * 
     * Simplexes with fewer than N-1 dimensions cannot be intersected, and the interval will
     * always be empty.
     */
    Rect<T,1> intersect(const Ray<T,N>& r) const {
        if (n == N + 1) {
            // where are barycentric coords equal to ray coords?
            // Σ w_i * p_i = o + s * v
            // Σ w_i       = 1
            // combined:
            // Σ (p_i, 1) * w_i = (o, 1) + s * (v, 0)
            // M * w            = B * (1, s)
            Vec<T,N+1> m[N+1];
            Vec<T,N+1> B[2] = {Vec<T,N+1>(r.origin, 1), Vec<T,N+1>(r.direction, 0)};
            for (index_t i = 0; i < N + 1; i++) {
                m[i] = Vec<T,N+1>(pts[i], 1);
            }
            // todo: handle degenerate simplexes
            // put the solutions (p_i, k_i) for `w_i = p_i + k_i * s` into B:
            if (linear_solve(m, 2, B)) {
                // intersect all the ranges of s for which 0 <= w_i <= 1.
                Rect<T,1> interval = Rect<T,1>::full;
                for (index_t i = 0; i < N + 1 and interval.hi > interval.lo; ++i) {
                    // w_i = p_i + k_i * s
                    // solve for w_i = 0 and w_i = 1:
                    T p_i = B[0][i];
                    T k_i = B[1][i];
                    if (k_i != 0) {
                        T s0 =     -p_i  / k_i;
                        T s1 = (1 - p_i) / k_i;
                        interval &= Rect<T,1>::from_corners(s0, s1);
                    } else if (p_i < 0 or p_i > 1) {
                        // if p_i is outside (0, 1) then no value of `s` can solve w_i in
                        // (0, 1), because with k = 0 there is no dependency on s!
                        return Rect<T,1>::empty;
                        // (otherwise, all values of `s` solve w_i in (0, 1) for this w_i;
                        // no change to the solution interval).
                    }
                }
                return interval;
            }
        } else if (n == N) {
            T s;
            if (trace_simplex<T,N>(pts, r, nullptr, &s)) {
                return Rect<T,1>(s, s);
            }
        }
        // no hit
        return Rect<T,1>();
    }
    
    /**
     * @brief Project the point to the simplex along its orthogonal complement
     * (if there is one) and test for containment.
     * 
     * @param p A point.
     * @param params Optional return variable; return the `n - 1` surface parameters
     * of this simplex representing the projection of `p` onto the space
     * spanned by this simplex.
     * @return `true` if `p` is on or inside this simplex after projection; 
     * `false` otherwise.
     */
    bool is_facing(const Vec<T,N>& p, Vec<T,N>* params=nullptr) const {
        switch (n) {
            case 0:
                return false;
            
            case 1: {
                return true;
            } break;
            
            case 2: {
                Vec<T,N> v = pts[1] - pts[0];
                T d = v.dot(p - pts[0]) / v.mag2();
                if (params) (*params)[0] = d;
                return d >= 0 and d <= 1;
            } break;
            
            default: {
                // simplex spans a subspace. orthogonally project p along its nullspace.
                // if the simplex spans the domain, then nullspace() early exits and
                // everything else works as normal.
                
                // nullspace comes first (cause we won't solve for its coords),
                // followed by the simplex basis.
                const index_t n_bases = n - 1;
                const index_t n_null  = N - n_bases;
                Vec<T,N> null_bases[N];
                Vec<T,N>* const bases = null_bases + n_null; 
                
                // compute simplex bases
                for (index_t i = 1; i < n; ++i) {
                    bases[i - 1] = pts[i] - pts[0];
                }
                
                // find the simplex's nullspace
                nullspace(bases, n_bases, null_bases);
                
                // write P as a linear combination of the simplex basis and 
                // its nullspace bases. (but don't bother actually solving for 
                // the nullpsace coords; we don't use them).
                Vec<T,N> x = p - pts[0];
                if (not linear_solve(null_bases, 1, &x, n_null)) return false;
                
                // check that the coordinates of `p` are inside the simplex
                T sum = 0;
                for (index_t i = n_null; i < n_bases; ++i) {
                    sum += x[i];
                    if (x[i] < 0 or sum > 1) return false;
                }
                
                // provide surface parameters to caller, e.g. so they can perform
                // the projection if they want to.
                if (params) {
                    for (index_t i = 0; i < n_bases; ++i) {
                        (*params)[i] = x[i + n_null];
                    }
                }
                
                return true;
            }
        }

        // statically unreachable, but some less-clever compilers may complain without:
        return false;
    }
    
    /**
     * @brief Compute the signed volume of the simplex.
     * 
     * If the simplex is not a full volume (i.e., the number of vertices is less than `N+1`),
     * then its volume is zero.
     */
    T measure_interior() const {
        if (n < N + 1) { return 0; }
        T k = 1;
        Vec<T,N> vs[N];
        for (index_t i = 1; i <= N; ++i) {
            vs[i - 1] = pts[i] - pts[0];
            k *= i;
        }
        T* m = vs[0].begin();
        if constexpr (N == 2) {
            return det2x2(m) / 2;
        } else if constexpr (N == 3) {
            return det3x3(m) / 6;
        } else {
            return det_destructive(m, N) / k;
        }
    }
    
    /**
     * @brief Return the positive volume of the simplex within the subspace that it spans.
     * 
     * For example, if the simplex has three points, it spans a plane and this function
     * returns the area within the plane. If it has four, it returns the volume of the
     * tetrahedron; or if it has two, the length of the line segment; and so on.
     * 
     * @return T 
     */
    T measure() const {
        // special cases (easy to compute):
        // line segment
        if (n == 2) return (pts[1] - pts[0]).mag();
        // 3D triangle
        if constexpr (N == 3) {
            if (n == 3) return (pts[1] - pts[0]).cross(pts[2] - pts[0]).mag() / 2;
        }
        // full volume
        if (n == N + 1) return measure_interior();
        // otherwise, construct the basis of the simplex.
        // we'll use it to compute the volume
        Vec<T,N> bases[N];
        int k = 1;
        for (index_t i = 1; i < n; ++i) {
            bases[i - 1] = pts[i] - pts[0];
            k *= i;
        }
        T* m = bases[0].begin();
        // volume is sqrt of the determinant of the gram matrix,
        // divided by n factorial. i suspect the gram matrix could probably
        // be computed more efficiently, because it's symmetric
        T buf[N * N];
        WrapperMatrix<T,0,0> b {buf, n, n};
        WrapperMatrix<T,0,N,MatrixLayout::ROW_MAJOR> m0 {m, n};
        WrapperMatrix<T,N,0,MatrixLayout::COL_MAJOR> m1 {m, n};
        // (n x n) = (n x N) * (N x n)
        mul(&b, m0, m1);
        return std::sqrt(std::abs(det_destructive(buf, n))) / k;
    }
    
    /**
     * @brief Simplex-point intersection test.
     * 
     * If this simplex has fewer than `N + 1` vertices, this simplex
     * is not a volume, and this function returns `false`.
     * 
     * If the simplex is degenerate (two or more coincident points),
     * then it is also considered empty and this function returns `false`.
     * 
     * @param p A point.
     * @return `true` if `p` is on or inside this simplex; `false` otherwise.
     */
    bool contains(const Vec<T,N>& p) const {
        if (n < N + 1) return false;
        // solve for barycentric coordinates x:
        // Σ (1, v[i]) * x[i] = (1, p)
        
        const index_t K = N + 1;
        T m[K * K];
        T x[K];
        for (index_t i = 0; i < K; ++i) {
            const T* v = pts[i].begin();
            m[i * K] = 1;
            std::copy(v, v + N, m + i * K + 1);
        }
        x[0] = 1;
        std::copy(p.begin(), p.end(), x + 1);
        
        // (no solution means a degenerate simplex, which has zero volume, so empty)
        if (not linear_solve<T,false>((T*) m, K, 1, (T*) x, 1)) return false;
        // skip solving for the sum of the weights         ⤴︎
        // we don't need it, and we know it's 1.
        
        T sum = 0;
        for (index_t i = 1; i < K; ++i) {
            sum += x[i];
            if (x[i] < 0 or sum > 1) return false;
        }
        return true;
    }
    
    
    /// Return the signed distance to the surface of the shape.
    T sdf(Vec<T,N> p, Vec<T,N>* normal=nullptr) const {
        detail::SimplexProjection<T,N> proj {
            *this,
            p,
            detail::ProjectionOp::PROJECT
        };
        if (normal) *normal = proj.normal_direction();
        return proj.distance() * (proj.result.contains ? -1 : 1);
    }
    
    Vec<T,N> normal(Vec<T,N> p) const {
        detail::SimplexProjection<T,N> proj {
            *this,
            p,
            detail::ProjectionOp::PROJECT
        };
        return proj.normal_direction();
    }
    
    Vec<T,N> project(Vec<T,N> p, Simplex<T,N>* face=nullptr) const {
        detail::SimplexProjection<T,N> proj {
            *this,
            p,
            detail::ProjectionOp::PROJECT
        };
        if (face) *face = proj.projected_face();
        return proj.projected_point();
    }
    
    /**
     * @brief Ensure `p` lies within the simplex by orthogonally projecting it to
     * the nearest point on the surface if it lies outside.
     *
     * The point is unchanged if `p` is already inside the simplex.
     * 
     * @param p The point to project onto the surface of this simplex.
     * @param onto Optional output parameter to receive the sub-simplex
     * onto which `p` was projected. It is permissible for `onto` to 
     * alias `this`.
     *
     * @return The location of `p`'s projection.
     */
    Vec<T,N> clip(const Vec<T,N>& p, Simplex<T,N>* onto=nullptr) const {
        detail::SimplexProjection<T,N> proj {
            *this,
            p,
            detail::ProjectionOp::CLIP
        };
        if (onto) *onto = proj.projected_face();
        return proj.projected_point();
    }
    
    void exclude(index_t i) {
        std::copy(pts + i + 1, pts + n, pts + i);
        n -= 1;
    }
    
    /**
     * @brief Create a new sub-simplex by excluding the `i`th vertex in this simplex.
     */
    Simplex<T,N> excluded(index_t i) const {
        Simplex<T,N> s;
        index_t c = 0;
        for (index_t j = 0; j < n; ++j) {
            if (i != j) s.pts[c++] = pts[j];
        }
        s.n = n - 1;
        return s;
    }
    
    Vec<T,N> convex_support(Vec<T,N> d) const {
        Vec<T,N> o = pts[0];
        T m = pts[0].dot(d);
        for (index_t i = 1; i < n; ++i) {
            T t = pts[i].dot(d);
            if (t > m) {
                m = t;
                o = pts[i];
            }
        }
        return o;
    }
    
    Rect<T,N> bounds() const {
        Rect<T,N> r = Rect<T,N>::empty;
        for (index_t i = 0; i < n; ++i) {
            r |= pts[i];
        }
        return r;
    }
    
    template <ConvexObject Shape>
    requires (Shape::N == N) and std::same_as<T, typename Shape::elem_t>
    bool intersects(const Shape& other) const {
        return geom::intersects(
            as_any_convex(*this),
            as_any_convex(other)
        );
    }

}; // class simplex

/// @addtogroup shape
/// @{

/// @brief Transform a simplex by a similarity transform.
/// @related Simplex
/// @related Similarity
template <typename T, index_t N>
Simplex<T,N> operator*(const Similarity<T,N>& xf, const Simplex<T,N>& s) {
    Simplex<T,N> s2 = s;
    for (index_t i = 0; i < s2.n; ++i) {
        s2.pts[i] = xf * s2.pts[i];
    }
    return s2;
}

/// @brief Inverse-transform a simplex by a similarity transform.
/// @related Simplex
/// @related Similarity
template <typename T, index_t N>
Simplex<T,N> operator/(const Simplex<T,N>& s, const Similarity<T,N>& xf) {
    Simplex<T,N> s2 = s;
    for (index_t i = 0; i < s2.n; ++i) {
        s2.pts[i] = s2.pts[i] / xf;
    }
    return s2;
}

/// @brief Transform a simplex by an isometry.
/// @related Simplex
/// @related Isometry
template <typename T, index_t N>
Simplex<T,N> operator*(const Isometry<T,N>& xf, const Simplex<T,N>& s) {
    Simplex<T,N> s2 = s;
    for (index_t i = 0; i < s2.n; ++i) {
        s2.pts[i] = xf * s2.pts[i];
    }
    return s2;
}

/// @brief Inverse-transform a simplex by an isometry.
/// @related Simplex
/// @related Isometry
template <typename T, index_t N>
Simplex<T,N> operator/(const Simplex<T,N>& s, const Isometry<T,N>& xf) {
    Simplex<T,N> s2 = s;
    for (index_t i = 0; i < s2.n; ++i) {
        s2.pts[i] = s2.pts[i] / xf;
    }
    return s2;
}

/**
 * Ray trace an N-1 dimensional simplex (e.g. triangle in 3D, tetrahedron in 4D, a line
 * in 2D etc.), returning the surface parameters of the hit point, given in coordinates
 * of the basis spanned by the edges radiating from `verts[0]`.
 * 
 * @param verts An array of the simplex's N vertices.
 * @param ray The ray to intersect with the simplex.
 * @param uv Buffer for the surface coordinates of the hitpoint. (return variable; may
 * be null).
 * @param s The ray parameter of the hit point. (return variable)
 * @return `true` if the ray hit the simplex.
 */
template <typename T, index_t N>
inline bool trace_simplex(const Vec<T,N> verts[N], const Ray<T,N>& ray, Vec<T,N-1>* uv, T* s) {
    /* Solve for u,v coordinates along the edges radiating from `verts[0]`:
    A + u(B-A) + v(C-A)      = o + sV
        u(B-A) + v(C-A) - sV = o - A
    Mx = o - A
    */
    
    // todo: cramer's rule might be used in 3d for faster linear solve. 
    //       (note: is unstable).
    
    // compute bases
    Vec<T,N> bases[N];
    for (index_t i = 1; i < N; ++i) {
        bases[i - 1] = verts[i] - verts[0];
    }
    bases[N - 1] = -ray.direction;
    
    // linear solve
    Vec<T,N> x = ray.origin - verts[0];
    index_t skip = uv ? 0 : N - 1;
    if (not linear_solve(bases, 1, &x, skip)) return false;
    
    // inside simplex?
    T sum = 0;
    for (index_t i = 0; i < N - 1; ++i) {
        if (x[i] < 0) return false;
        sum += x[i];
    }
    if (sum > 1) return false;
    
    // output result
    *s  = x[N];
    if (uv) *uv = x.template resized<N-1>();
    return true;
}

template <typename T>
using Triangle = Simplex<T,2>;

template <typename T>
using Tetrahedron = Simplex<T,3>;

/// @} // group shape

template <typename T, index_t N, typename H>
struct Digest<Simplex<T,N>, H> {
    H operator()(const Simplex<T,N>& s) const {
        H nonce = geom::truncated_constant<H>(0x8159bb983f797983, 0xe8904e94430f7555);
        return geom::hash_array<Vec<T,N>,H>(nonce, s.pts, s.n);
    }
};

#ifdef GEOMC_USE_STREAMS

template <typename T, index_t N>
std::ostream& operator<<(std::ostream& os, const Simplex<T,N>& s) {
    os << "Simplex(";
    for (index_t i = 0; i < s.n; ++i) {
        if (i > 0) os << ", ";
        os << s.pts[i];
    }
    os << ")";
    return os;
}

#endif

} // namespace geom


template <typename T, index_t N>
struct std::hash<geom::Simplex<T,N>> {
    size_t operator()(const geom::Simplex<T,N> &s) const {
        return geom::hash<geom::Simplex<T,N>, size_t>(s);
    }
};
