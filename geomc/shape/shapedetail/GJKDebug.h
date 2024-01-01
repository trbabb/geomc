#pragma once

#include <geomc/shape/Simplex.h>

namespace geom {

/* debug file format:
 * 
 * convex hull points of shape A
 * convex hull points of shape B
 * [
 *   simplex before projection
 *   simplex after projection
 *   pt on minkowski difference
 *   search direction from simplex toward origin
 *   ---
 * ] x number of steps in the search
 * ========
 */


template <typename T, index_t N>
void emit_simplex(FILE* f, const Simplex<T,N>& s) {
    fprintf(f, "simplex %ldD %ld | ", N, s.n);
    for (index_t i = 0; i < s.n; ++i) {
        for (int j = 0; j < N; ++j) {
            fprintf(f, "  %f", s.pts[i][j]);
        }
        fprintf(f, ", ");
    }
    fprintf(f, "\n");
}

template <typename T, index_t N>
void emit_vec(FILE* f, const char* label, const Vec<T,N>& v) {
    fprintf(f, "%s %ldD | ", label, N);
    for (int j = 0; j < N; ++j) {
        fprintf(f, "  %f", v[j]);
    }
    fprintf(f, "\n");
}

template <typename T, index_t N>
void emit_splex_stage(
        FILE* f,
        const Simplex<T,N>& cur,
        const Simplex<T,N>& next,
        const Vec<T,N>& a,
        const Vec<T,N>& d)
{
    emit_simplex(f, cur);
    emit_simplex(f, next);
    emit_vec(f, "A", a); // we don't actually need to emit this, because it'll be the last vtx in `cur`
    emit_vec(f, "p", d);
    fprintf(f, "---\n");
}

template <typename T, index_t N>
void emit_hull(FILE* f, const Vec<T,N>* pts, index_t n) {
    fprintf(f, "hull %ldD %ld | ", N, n);
    for (index_t i = 0; i < n; ++i) {
        for (int j = 0; j < N; ++j) {
            fprintf(f, "  %f", pts[i][j]);
        }
        fprintf(f, ", ");
    }
    fprintf(f, "\n");
}


} // namespace geom
