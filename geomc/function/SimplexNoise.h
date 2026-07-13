#pragma once
/*
 * SimplexNoise.h
 *
 *  Created on: Jul 12, 2026
 *      Author: tbabb
 */

#include <array>
#include <utility>

#include <geomc/function/Utils.h>
#include <geomc/random/SampleVector.h>
#include <geomc/random/DenseDistribution.h>
#include <geomc/linalg/Vec.h>
#include <geomc/linalg/Matrix.h>

namespace geom {

/**
 * @ingroup function
 * @brief Real-valued smooth simplex noise over `N` dimensions.
 *
 * Like PerlinNoise, this produces a continuous, band-limited scalar field with
 * random unit gradients placed on a lattice. Unlike PerlinNoise — which
 * interpolates the _O(2<sup>N</sup>)_ corners of a hypercube — simplex noise
 * sums the contributions of the _N + 1_ corners of the single lattice simplex
 * containing the sample point. This makes it _O(N<sup>2</sup>)_ per sample
 * rather than _O(2<sup>N</sup>)_, which is dramatically cheaper in high
 * dimensions (the two are roughly equal cost near `N = 4`).
 *
 * The gradients are drawn uniformly from the unit sphere — the same
 * high-quality gradients used by PerlinNoise — rather than the small fixed
 * gradient set used by many simplex noise implementations.
 *
 * The field is C<sup>1</sup> continuous: the interpolation kernel and its
 * derivative both vanish exactly where a corner leaves the containing simplex,
 * so the analytic gradient returned by `gradient()` is continuous everywhere.
 */
template <typename T, index_t N>
class SimplexNoise : public Dimensional<T,N> {
    static constexpr size_t N_GRADIENTS = 0x100;

public:

    using gridtype  = PointType<index_t,N>;
    using pointtype = PointType<T,N>;
    using grid_t    = typename gridtype::point_t;
    using typename Dimensional<T,N>::point_t;

private:

    std::shared_ptr<point_t[]> _gradients;
    T _scale; // amplitude normalization, so the output lies in roughly [-1, 1]

public:

    /**
     * Construct a new simplex noise object with a `std::random_device` as a
     * source of random bits.
     */
    SimplexNoise():SimplexNoise(std::random_device()) {}

    /**
     * Construct a new simplex noise object with the supplied random number generator.
     * @param rng A source of random bits.
     */
    template <typename Generator>
    SimplexNoise(Generator& rng): _gradients(new point_t[N_GRADIENTS]) {
        // gradients are random unit vectors (for N == 1, random signs).
        for (index_t i = 0; i < (index_t)N_GRADIENTS; ++i) {
            _gradients[i] = random_unit<T,N>(rng);
        }
        _scale = amplitude_scale();
    }

    T operator()(point_t pt) const {
        return eval(pt);
    }

    /**
     * Evaluate the noise at `pt`.
     */
    T eval(point_t pt) const {
        return _scale * raw_eval(_gradients.get(), pt);
    }

    /**
     * Evaluate the gradient of the noise function at `pt`.
     *
     * @param pt Location at which to sample the noise function.
     * @return A pair of (`noise(x)`, `gradient(noise(x))`).
     */
    std::pair<T, point_t> gradient(point_t pt) const {
        const point_t* grads = _gradients.get();
        T       val  = 0;
        point_t grad = {};
        for_each_corner(skew(pt), [&](const grid_t& corner) {
            point_t g = grads[hash_index(corner)];
            point_t d = pt - unskew(corner);
            T w = kernel_radius2() - pointtype::mag2(d);
            if (w > 0) {
                // contribution is  w^4 * (g . d),  with  w = r^2 - |d|^2.
                // d/dx of that is  -8 w^3 (g . d) d  +  w^4 g.
                T gd = pointtype::dot(g, d);
                T w2 = w * w;
                T w3 = w2 * w;
                val  += w2 * w2 * gd;
                grad += (-8 * w3 * gd) * d + (w2 * w2) * g;
            }
        });
        return std::pair<T, point_t>(_scale * val, _scale * grad);
    }

    /**
     * @brief The matrix that maps the integer lattice to the skewed simplex
     * lattice in sample space (an equilateral parallelepiped).
     *
     * The columns are the edge vectors of the parallelepiped; they have equal
     * length and equal pairwise angles. This is the inverse of `skew_matrix()`.
     */
    static const SimpleMatrix<T,N,N>& unskew_matrix() {
        // M = I - G * J, where J is the all-ones matrix. Applying it to a lattice
        // point is done directly as a rank-1 update in `unskew()`; this matrix is
        // provided for completeness and for callers who want the transform.
        static const SimpleMatrix<T,N,N> m = build_skew_matrix(-unskew_factor());
        return m;
    }

    /**
     * @brief The matrix that maps sample space to the integer lattice, skewing
     * the simplex lattice onto the integer grid so the containing cell can be
     * found by flooring. This is the inverse of `unskew_matrix()`.
     */
    static const SimpleMatrix<T,N,N>& skew_matrix() {
        // M^-1 = I + F * J.
        static const SimpleMatrix<T,N,N> m = build_skew_matrix(skew_factor());
        return m;
    }

protected:

    // raw (un-normalized) noise value, using an explicit gradient table.
    // shared by eval() and the amplitude-normalization estimate.
    static T raw_eval(const point_t* grads, point_t pt) {
        T sum = 0;
        for_each_corner(skew(pt), [&](const grid_t& corner) {
            point_t d = pt - unskew(corner);
            T w = kernel_radius2() - pointtype::mag2(d);
            if (w > 0) {
                // kernel is w^4; it falls smoothly to zero in both value and
                // slope at the simplex boundary, where this corner leaves the
                // containing simplex. that is what makes the field C^1.
                T w2 = w * w;
                sum += w2 * w2 * pointtype::dot(grads[hash_index(corner)], d);
            }
        });
        return sum;
    }

    // hash an integer lattice vertex to a gradient table index.
    static index_t hash_index(const grid_t& corner) {
        uint64_t idx = 0;
        auto m = gridtype::iterator(corner);
        for (index_t i = 0; i < N; ++i) {
            // iterated linear congruential scramble, using Knuth's constants:
            uint64_t k = static_cast<uint64_t>(m[i]);
            idx = 6364136223846793005ULL * (k + idx) + 1442695040888963407ULL;
        }
        return idx & (N_GRADIENTS - 1);
    }

    // skew constants. world -> lattice via (I + F J); lattice -> world via (I - G J),
    // where J is the all-ones matrix. These are the standard simplex-noise skew
    // factors, chosen so the unit cube maps to an equilateral parallelepiped.
    static T skew_factor() {
        return (std::sqrt((T)N + 1) - 1) / (T)N;
    }

    static T unskew_factor() {
        return (1 - 1 / std::sqrt((T)N + 1)) / (T)N;
    }

    // apply (I + F J) to a sample point: x_i + F * sum(x). O(N) rank-1 form of
    // (skew_matrix() * pt).
    static point_t skew(point_t pt) {
        return pt + skew_factor() * component_sum(pt);
    }

    // apply (I - G J) to a lattice point: c_i - G * sum(c). O(N) rank-1 form of
    // (unskew_matrix() * corner).
    static point_t unskew(const grid_t& corner) {
        point_t p;
        T s = 0;
        auto m  = gridtype::iterator(corner);
        auto pi = pointtype::iterator(p);
        for (index_t i = 0; i < N; ++i) {
            T c   = (T)m[i];
            pi[i] = c;
            s    += c;
        }
        return p - unskew_factor() * s;
    }

    static T component_sum(point_t pt) {
        T s = 0;
        auto m = pointtype::iterator(pt);
        for (index_t i = 0; i < N; ++i) s += m[i];
        return s;
    }

    /**
     * @brief Squared radius of the interpolation kernel.
     *
     * Equal to the minimum simplex height, squared: the largest ball centered on
     * a lattice vertex that fits inside the union of simplices touching it. With
     * this radius, a corner's contribution vanishes (in both value and gradient)
     * exactly as the sample point crosses out of the containing simplex, which is
     * what makes the field C<sup>1</sup>.
     */
    static T kernel_radius2() {
        static const T r2 = compute_kernel_radius2();
        return r2;
    }

    // Freudenthal / Kuhn triangulation: enumerate the N + 1 lattice vertices of
    // the simplex containing the (skewed) point `u`, calling `fn(grid_t)` for each.
    // The simplex is the path from the base corner to the opposite corner that
    // steps along the axes in order of decreasing fractional coordinate.
    template <typename Fn>
    static void for_each_corner(point_t u, Fn&& fn) {
        grid_t  base;
        point_t frac;
        auto ub = pointtype::iterator(u);
        auto bi = gridtype::iterator(base);
        auto fi = pointtype::iterator(frac);
        for (index_t i = 0; i < N; ++i) {
            T f   = std::floor(ub[i]);
            bi[i] = (index_t)f;
            fi[i] = ub[i] - f;
        }
        // sort axes by decreasing fractional coordinate (insertion sort; N is small).
        std::array<index_t,N> order;
        for (index_t i = 0; i < N; ++i) order[i] = i;
        for (index_t i = 1; i < N; ++i) {
            index_t a = order[i];
            index_t j = i;
            while (j > 0 and fi[order[j - 1]] < fi[a]) {
                order[j] = order[j - 1];
                j -= 1;
            }
            order[j] = a;
        }
        // walk the simplex corners, stepping one axis at a time.
        grid_t corner = base;
        auto ci = gridtype::iterator(corner);
        for (index_t k = 0; k <= N; ++k) {
            fn(corner);
            if (k < N) ci[order[k]] += 1;
        }
    }

    // build I + factor * J as an explicit matrix.
    static SimpleMatrix<T,N,N> build_skew_matrix(T factor) {
        SimpleMatrix<T,N,N> m; // SimpleMatrix default-constructs to the identity
        for (index_t r = 0; r < N; ++r) {
            for (index_t c = 0; c < N; ++c) {
                m(r,c) += factor;
            }
        }
        return m;
    }

    // kernel radius^2 = min over the reference simplex's vertices of the squared
    // distance from that vertex to the affine hull of the opposite facet.
    static T compute_kernel_radius2() {
        // reference Kuhn simplex in lattice space: 0, e0, e0+e1, ..., (1,1,...,1),
        // mapped into sample space.
        std::array<point_t, N + 1> w;
        grid_t corner = {};
        auto ci = gridtype::iterator(corner);
        for (index_t k = 0; k <= N; ++k) {
            w[k] = unskew(corner);
            if (k < N) ci[k] += 1;
        }
        T best = std::numeric_limits<T>::max();
        for (index_t k = 0; k <= N; ++k) {
            // height of w[k] above the affine hull of the other N vertices:
            // Gram-Schmidt the facet edges, then take the residual of
            // (w[k] - w[a]) orthogonal to the facet.
            index_t a = (k == 0) ? 1 : 0;
            std::array<point_t, N> facet_basis;
            index_t n_basis = 0;
            for (index_t j = 0; j <= N; ++j) {
                if (j == k or j == a) continue;
                point_t v = w[j] - w[a];
                for (index_t b = 0; b < n_basis; ++b) {
                    v = v - pointtype::dot(v, facet_basis[b]) * facet_basis[b];
                }
                T len = std::sqrt(pointtype::mag2(v));
                if (len > 0) facet_basis[n_basis++] = v / len;
            }
            point_t r = w[k] - w[a];
            for (index_t b = 0; b < n_basis; ++b) {
                r = r - pointtype::dot(r, facet_basis[b]) * facet_basis[b];
            }
            best = std::min(best, pointtype::mag2(r));
        }
        return best;
    }

    /**
     * @brief Amplitude normalization, computed once per `<T,N>`.
     *
     * Uses a fixed reference gradient table so all instances share the same
     * scale; the maximum amplitude is essentially independent of the specific
     * gradients, so a shared constant is appropriate.
     */
    static T amplitude_scale() {
        static const T s = compute_amplitude_scale();
        return s;
    }

    static T compute_amplitude_scale() {
        DefaultLCG rng(0x9E3779B97F4A7C15ULL);
        std::shared_ptr<point_t[]> grads(new point_t[N_GRADIENTS]);
        for (index_t i = 0; i < (index_t)N_GRADIENTS; ++i) {
            grads[i] = random_unit<T,N>(rng);
        }
        DenseUniformDistribution<T> u(-64, 64);
        T mx = 0;
        for (index_t s = 0; s < 8192; ++s) {
            point_t pt;
            auto pi = pointtype::iterator(pt);
            for (index_t i = 0; i < N; ++i) pi[i] = u(rng);
            mx = std::max(mx, std::abs(raw_eval(grads.get(), pt)));
        }
        return mx > 0 ? (T)1 / mx : (T)1;
    }
};

} // namespace geom
