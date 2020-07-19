#include <utility>

#include <geomc/shape/Rect.h>

namespace geom {

namespace detail {

/*****************************************
 * KDTree helpers                        *
 *****************************************/


// basic shape
// (works on rects)
template <typename T, index_t N, typename D>
struct ShapeIndexHelper {
    
    // can be anything unionable with a Rect
    // (i.e. either a rect or a vec)
    typedef Rect<T,N> bound_t;
    
    static inline Vec<T,N> getPoint(const D& obj) {
        return obj.center();
    }
    
    static inline T dist2(const D& obj, const Vec<T,N>& p) {
        return obj.dist2(p);
    }
    
    static inline bound_t bounds(const D& obj) {
        return obj.bounds();
    }
    
};


// pointers to indexable shapes
template <typename T, index_t N, typename D>
struct ShapeIndexHelper< T, N, D* > {
    
    typedef typename ShapeIndexHelper<T,N,D>::bound_t bound_t;
    
    static inline Vec<T,N> getPoint(const D* obj) {
        return ShapeIndexHelper<T,N,D>::getPoint(*obj);
    }
    
    static inline T dist2(const D* obj, const Vec<T,N>& p) {
        return ShapeIndexHelper<T,N,D>::getDist2(*obj, p);
    }
    
    static inline bound_t bounds(const D* obj) {
        return ShapeIndexHelper<T,N,D>::bounds(*obj);
    }
    
};


// key / value pairs
// (delegate to the shape helper of the key)
template <typename T, index_t N, typename K, typename V>
struct ShapeIndexHelper< T, N, std::pair<K, V> > {
    
    typedef typename ShapeIndexHelper<T,N,K>::bound_t bound_t;
    
    static inline Vec<T,N> getPoint(const std::pair<K,V>& obj) {
        return ShapeIndexHelper<T,N,K>::getPoint(obj.first);
    }
    
    static inline T dist2(const std::pair<K,V>& obj, const Vec<T,N>& p) {
        return ShapeIndexHelper<T,N,K>::getDist2(obj.first);
    }
    
    static inline bound_t bounds(const std::pair<K,V>& obj) {
        return ShapeIndexHelper<T,N,K>::bounds(obj.first);
    }
    
};


// bare vector 
template <typename T, index_t N>
struct ShapeIndexHelper< T, N, Vec<T,N> > {
    
    typedef Vec<T,N> bound_t;
    
    static inline Vec<T,N> getPoint(const Vec<T,N>& v) {
        return v;
    }
    
    static inline T dist2(const Vec<T,N>& v, const Vec<T,N>& p) {
        return v.dist2(p);
    }
    
    static inline Vec<T,N> bounds(const Vec<T,N>& v) {
        return v;
    }
    
};


// bounded shapes
template <typename T, index_t N>
struct ShapeIndexHelper< T, N, Bounded<T,N> > {
    
    typedef Rect<T,N> bound_t;
    
    inline Vec<T,N> getPoint(const Bounded<T,N>& r) {
        return r.bounds().center();
    }
    
    static inline T dist2(const Bounded<T,N>& r, const Vec<T,N>& p) {
        return r.bounds().clamp(p).dist2(p);
    }
    
    static inline bound_t bounds(const Bounded<T,N>& r) {
        return r.bounds();
    }
    
};


// returns an upper bound on the squared distance to the nearest object
// to `p` in mimimum bounding box `r`.
template <typename T, index_t N>
T worst_case_nearest2(const Vec<T,N>& p, const Rect<T,N>& r) {
    T result = std::numeric_limits<T>::max();
    
    // each face of the box must have an object touching it, if it is a minimal box.
    // each face has a farthest corner where such an object could be hiding,
    // and we want the nearest of these.
    
    Vec<T,N> far_extremes;
    Vec<T,N> near_extremes;
    
    for (index_t axis = 0; axis < N; axis++) {
        T lodist = std::abs(p[axis] - r.lo[axis]);
        T hidist = std::abs(p[axis] - r.hi[axis]);
        
        far_extremes[axis]  = lodist > hidist ? r.lo[axis] : r.hi[axis];
        near_extremes[axis] = lodist > hidist ? r.hi[axis] : r.lo[axis];
    }
    
    for (index_t axis = 0; axis < N; axis++) {
        Vec<T,N> face_extreme = far_extremes;
        // along this axis, we only need the nearer face:
        face_extreme[axis] = near_extremes[axis];
        result = std::min(result, (face_extreme - p).mag2());
    }
    return result;
}
    
} // namespace detail

} // namespace geom