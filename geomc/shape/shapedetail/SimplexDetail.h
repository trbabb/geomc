#ifndef SIMPLEX_DETAIL_H
#define SIMPLEX_DETAIL_H

#include <geomc/linalg/Vec.h>
#include <geomc/shape/ShapeTypes.h>

namespace geom {
namespace detail {


// todo: I think we need to try again and use the Intersect formulation as a template,
//   but make it recursive (or iterative unrolled).

/**
 * Keep a view into a Simplex which consists of a subset of its vertices.
 * 
 * A SubSimplex steals a pointer into the Simplex's buffer; therefore the donor
 * Simplex *must not* be freed if any live SubSimplex refers to it. 
 */
template <typename T, index_t N>
struct SubSimplex {
    // todo: can all this be done with a bit vector?
    //  just need an iterator over the indices of set bits
    const Vec<T,N>* pts;
    index_t idxs[N];
    index_t n;
    
    SubSimplex(const Vec<T,N>* pts, index_t n):
        pts(pts),
        n(n) {}
    
    SubSimplex(const Simplex<T,N>* s):
            pts(s->pts),
            n(s->n)
    {
        for (index_t i = 0; i < n; ++i) {
            idxs[i] = i;
        }
    }
    
    SubSimplex<T,N> exclude(index_t i) {
        SubSimplex<T,N> other(pts, n - 1);
        for (index_t j = i; j < n - 1; ++j) {
            index_t src   = (j >= i) ? j + 1 : j;
            other.idxs[j] = idxs[src];
        }
        return other;
    }
    
    const Vec<T,N>& operator[](index_t i) const {
        return pts[idxs[i]];
    }
    
    Simplex<T,N> to_simplex() const {
        Simplex<T,N> s;
        for (index_t i = 0; i < n; ++i) {
            s.pts[i] = (*this)[i];
        }
        s.n = n;
        return s;
    }
    
};


// works, but marginally slower than existing Simplex implementation :(
//   - simpler, but requires a matrix solve of larger size, so probably not the right
//     formulation.
//   - also: a previous checkin has [vtxs] and [normals] flipped. it's about the same
//     speed but it's simpler.
template <typename T, index_t N, bool test_root_vtx>
index_t find_nearest_face_barycentric(
        index_t     n,          // number of pts in the simplex
        Vec<T,N+1>  basis[N+1], // homogeneous vectors, laid out like: [normals] [vtxs]
        Vec<T,N+1>  P,          // homogeneous point to test for containment
        Vec<T,N+1>* solution,   // output variable for the solution
        Vec<T,N+1>  v_inward,   // a vector pointing "inward" from the face
        bool        is_root=true)
{
    const index_t n_null = N - n + 1;
    if (n == 1) {
        if (v_inward.dot(P - basis[N]) > 0) {
            // P is on the "interior" side of this point
            return 0;
        }
        // projection to a point is trivial.
        *solution      = Vec<T,N+1>();
        (*solution)[N] = 1;
        return n;
    } else if (n == 2) {
        
    } else {
        T sign = 1;
        // the root simplex has had its nullspace solved already
        if (not is_root) {
            Vec<T,N+1>& normal = basis[0] = orthogonal(basis + 1);
            // make the normal into a direction vector.
            // algorithmically this is okay because we know the other points are
            // on the projective plane. (Less sure if this numerically okay).
            normal[N] = 0;
            // if the calculated normal faces toward the interior of this simplex, 
            // perform a sign flip.
            sign = (normal.dot(v_inward) > 0) ? -1 : 1;
        }
        
        // write P as a weighted combination of each vertex, plus possibly
        // some vector(s) orthogonal to the simplex. (we must make a copy
        // of the basis because linear_solve() screws with the contents).
        Vec<T,N+1> buffer[N + 1];
        std::copy(basis, basis + N + 1, buffer);
        bool ok = linear_solve(buffer, solution, P, n_null - 1);
        Vec<T,N+1>& x = *solution;
        
        // all of the barycentric coords must lie inside the simplex
        for (index_t i = n_null; ok and i <= N; ++i) {
            if (x[i] < 0 or x[i] > 1) {
                ok = false;
            }
        }
        
        if (ok) {
            // P projects to the interior of this face.
            if (n <= N and sign * x[0] < 0 and not is_root) {
                // P is on the "backside" of this face; this face does not own it.
                return 0;
            } else {
                // P is on the "frontside" of this face, and therefore projects to it.
                return n;
            }
        }
        
        // P does not project into the interior of this simplex, but this simplex
        // faces it, so it must belong to one of the sub-faces.
        // try excluding verts from this simplex until we find it.
        
        // basis[0] will be clobbered by the sub-simplex's normal, so
        // use it to replace the first excluded vertex. if it's a normal,
        // it becomes part of the nullspace. if not, it's ignored.
        Vec<T,N+1> prev_excluded = basis[0];
        const index_t skip = test_root_vtx ? 0 : 1;
        for (index_t i = n_null; i <= N - skip; ++i) {
            std::swap(basis[i], prev_excluded);
            Vec<T,N+1> inward_dir = basis[N] - prev_excluded;
            index_t sub_n = find_nearest_face_barycentric<T,N,test_root_vtx>(
                    n - 1, basis, P, solution, inward_dir, false);
            if (sub_n > 0) {
                // P projects to this sub-face; we are done.
                return sub_n;
            }
        }
        
        // P did not belong to us, somehow. probably a degeneracy
        // xxx todo: just return self
        return n;
    }
}


// xxx: switched off
template <typename T, index_t N, bool test_root_vtx>
bool find_nearest_face(
        SubSimplex<T,N> s,            // sub-simplex to test for projection containment
        Vec<T,N>        to_p,         // vec from a pt on the simplex to `p`
        Vec<T,N>        all_bases[N], // [spanning basis] [null basis]
        Vec<T,N>        outward_dir,  // a vector pointing "away" from the sub-simplex
        bool            is_root,      // is the sub-simplex the root simplex?
        Simplex<T,N>*   out)          // output simplex
{
    const index_t n_bases = s.n - 1;
    if (s.n == 1) {
        if (to_p.dot(outward_dir) < 0) return false;
        // projection onto a single point is trivial.
        // note: we leave `all_bases` in a garbage state,
        // but we don't use it, since we don't need it!
        *out = s.to_simplex();
        return true;
    } else {
        Vec<T,N> normal;
        if (not is_root) {
            // if this simplex is a face, is P inside or outside of it?
            // (the root simplex does not have an "inside" or an "outside")
            normal = orthogonal(all_bases);
            if ((to_p.dot(normal) < 0) == (normal.dot(outward_dir)) > 0) {
                // `p` falls "inward" of this face. `p` therefore does not project
                // to this simplex; instead it projects to some sibling sub-simplex.
                return false;
            }
        } else {
            // if the normal exists, it's been computed already. otherwise
            // the simplex spans the space, and the value of the "normal"
            // is irrelevant. We choose the final basis vector because orthogonal()
            // only examines (N - 1) elements of its input array, so this is 
            // the basis that would have been ignored:
            normal = all_bases[N - 1];
        }
        // test if P projects to any of our sub-simplexes.
        // try excluding each vertex, one at a time:
        Vec<T,N> prev_excluded = normal;
        index_t  prev_excl_idx = s.idxs[n_bases];
        s.n                   -= 1;
        for (index_t i = n_bases - 1; i >= 0; --i) {
            // put the previously-excluded basis back, and exclude a different one.
            // (the first basis we exclude will be replaced with our fresh null basis,
            // and remain there).
            std::swap(all_bases[i], prev_excluded);
            std::swap(s.idxs[i],    prev_excl_idx);
            if (find_nearest_face<T,N,test_root_vtx>(
                    s, to_p, all_bases, -prev_excluded, false, out))
            {
                // some sub-simplex of us both faces `p` and contains
                // its projection. `all_bases` and `out` now reflect that sub-simplex.
                return true;
            }
        }
        
        s.idxs[n_bases - 1] = prev_excl_idx;
        
        if (test_root_vtx) {
            // try excluding the basis-defining vertex. this means re-computing all the
            // basis vectors in terms of a new "root" vertex.
            for (index_t i = 0; i < n_bases - 1; ++i) {
                all_bases[i] = s[i] - s[n_bases - 1];
            }
            // a vector pointing from the new basis vertex to the old one:
            Vec<T,N> change_of_base = s[n_bases] - s[n_bases - 1];
            to_p += change_of_base;
            // test that sub-simplex.
            if (find_nearest_face<T,N,test_root_vtx>(
                    s, to_p, all_bases, -change_of_base, false, out))
            {
                return true;
            }
            // replace the normal with the remaining spanning basis.
            all_bases[n_bases - 1] = change_of_base;
        } else {
            // we left the normal in place of the last spanning basis.
            // put the spanning basis back.
            all_bases[n_bases - 1] = prev_excluded;
        }
        
        // complete the nullspace
        if (n_bases < N) {
            all_bases[N - 1] = normal;
        }
        
        // no sub-simplex claimed P was "beyond" it and took the projection.
        // therefore, this is the simplex onto which `p` projects.
        s.n += 1;
        *out = s.to_simplex();
        return true;
    }
}


} // namespace detail
} // namespace geom

#endif // define SIMPLEX_DETAIL_H
