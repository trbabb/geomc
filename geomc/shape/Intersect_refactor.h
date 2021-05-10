namespace geom {


// todo: we don't need to try excluding the vertex that was just added; 
//       that would revert to a simplex that was already tried.
//       if none of the remaining sub-simplices "contain" the origin,
//       then reject the intersection test.
    
// todo: does the algorithm break if we check the "back" face?
//       (It currently checks the back face).
//       (I think "no", but really not sure).

// todo: would be better to compute bases once, and then manipulate
//       ptrs to them. as dimensions get higher, this gets more important.
//       > should implement a helper class and delegate to it from Simplex<T,N>.
//       > we also currently have to keep re-computing the bases because we might exclude
//         the simplex origin. In GJK, we can construct the simplex origin
//         to be the last-added vertex (which cannot be excluded), allowing us
//         to use the same basis vectors for the entire calculation.
    
// todo: see Orthogonal.h— we compute L in many of our linear solves, but
//       we don't need it for nullspaces. that can save us time.

template <typename T, index_t N>
bool gjk_intersect(
    const Convex<T,N>& shape_a, 
    const Convex<T,N>& shape_b,
    Vec<T,N>* overlap_axis,
    Simplex<T,N>* simplex)
{
    Vec<T,N> d; // test separation axis
    Vec<T,N> a; // latest pt on minkowski difference
    Simplex<T,N>& s = *simplex;
    s.n = 0;
    
    // choose an arbitrary initial direction.
    if (overlap_axis and not overlap_axis.isZero()) {
        d = *overlap_axis;
    } else {
        d[0] = 1;
    }
    
    while (true) {
        a = shape_a.convex_support( d) - 
            shape_b.convex_support(-d);
        s |= a;
        
        if (a.dot(d) < 0) {
            if (overlap_axis) *overlap_axis = -a.projectOn(d);
            return false;
        }
        
        // find the sub-simplex facing the origin, and project the origin to it.
        // this vector is normal to the new simplex and connects the simplex's face to the origin.
        d = -s.project(Vec<T,N>(), &s);
        
        if (s.n == N + 1) {
            // full containment; no vertices were discarded by the projection.
            // the origin is inside the simplex.
            if (overlap_axis) *overlap_axis = d;
            return true;
        }
    }
}


} // end namespace geom