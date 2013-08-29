/*
 * AffineTransform.h
 *
 *  Created on: Oct 11, 2009
 *      Author: tbabb
 */

//TODO: This should not accept DYNAMIC_DIM.

#ifndef AFFINETRANSFORM_H_
#define AFFINETRANSFORM_H_

#include <algorithm>

#include <geomc/linalg/LinalgTypes.h>
#include <geomc/linalg/Vec.h>
#include <geomc/linalg/Quaternion.h>
#include <geomc/linalg/Ray.h>
#include <geomc/linalg/Matrix.h>

#ifdef GEOMC_LINALG_USE_STREAMS
    #include <iostream>
#endif

namespace geom {
    
/** @ingroup linalg
 *  @brief Affine transformation class.
 *
 * Vectors are assumed to be columns; therefore a transformation is ordered like:
 * 
 *     T * v
 * 
 * and so:
 * 
 *     A * B * v
 * 
 * applies B to v, then A to the result.
 * 
 * The transform which "undoes" a transform can be obtained with `xf.inverse()`. 
 * Because an AffineTransform internally stores both itself and its inverse, this 
 * construction is fairly inexpensive, as are all inverse-application functions.
 */
    
    template <typename T, index_t N> class AffineTransform {
    public:
        /// Matrix representing this transformation.
        SimpleMatrix<T,N+1,N+1> mat;
        /// Matrix representing the inverse of this transformation.
        SimpleMatrix<T,N+1,N+1> inv;

        /// Construct a new identity transform.
        AffineTransform() {}
        
        /*******************************
         * Operators                   *
         *******************************/
        
        /**
         * Transformation of a ray.
         */
        friend Ray<T,N> operator*(const AffineTransform<T,N> &at, Ray<T,N> r) {
            return Ray<T,N>(at.apply(r.origin), at.applyVector(r.direction));
        }
        
        /**
         * Inverse transformation of a ray (`xf`<sup>`-1`</sup>` * ray`)
         */
        friend Ray<T,N> operator/(Ray<T,N> r, const AffineTransform<T,N> &at) {
            return Ray<T,N>(at.applyInverse(r.origin), at.applyInverseVector(r.direction));
        }
        
        /**
         * Transformation of a point.
         */
        friend Vec<T,N> operator*(const AffineTransform<T,N> &at, Vec<T,N> p) {
            return at.apply(p);
        }
        
        /**
         * Inverse transformation of a point (`xf`<sup>`-1`</sup>` * pt`)
         */
        friend Vec<T,N> operator/(Vec<T,N> p, const AffineTransform<T,N> &at) {
            return at.applyInverse(p);
        }
        
        /**
         * Concatenation of transforms. 
         * @return A transformation representing an application of `at2` followed by
         * `at1`.
         */
        friend AffineTransform<T,N> operator*(const AffineTransform<T,N> &at1, const AffineTransform<T,N> &at2) {
            return at2.apply(at1);
        }
        
        /**
         * Concatenation of transforms.
         * 
         * Assign a transform representing an application of `this` followed by `at`.
         */
        AffineTransform<T,N>& operator*=(const AffineTransform<T,N> &at) {
            mat = at.mat * mat;
            inv = inv * at.inv;
            return (*this);
        }
        
        #ifdef GEOMC_LINALG_USE_STREAMS
        friend std::ostream &operator<<(std::ostream &stream, const AffineTransform<T,N> at) {
            stream << at.mat;
            return stream;
        }
        #endif
        
        /*******************************
         * Methods                     *
         *******************************/
        
        /**
         * Transformation of a point.
         */
        const Vec<T,N> apply(const Vec<T,N> &p) const {
            Vec<T,N+1> p_hom(p,1);
            p_hom = mat * p_hom;
            return p_hom.template resized<N>(); //c++ is awful
        }
        
        /**
         * Transformation of a direction vector; ignores any translation.
         */
        const Vec<T,N> applyVector(const Vec<T,N> &v) const {
            Vec<T,N> o;
            for (index_t r = 0; r < N; r++) {
                for (index_t c = 0; c < N; c++) {
                    o[r] += mat.get(r,c) * v[c];
                }
            }
            return o;
        }
        
        /**
         * Transformation of a normal. Preserves surface direction
         * of geometry transformed by `this`. 
         */
        const Vec<T,N> applyNormal(const Vec<T,N> &n) const {
            // normal matrix = txpose of inverse.
            Vec<T,N> o;
            for (index_t r = 0; r < N; r++) {
                for (index_t c = 0; c < N; c++) {
                    o[r] += inv.get(c,r) * n[c];
                }
            }
            return o;
        }
        
        /**
         * Inverse transformation of a point.
         */
        const Vec<T,N> applyInverse(const Vec<T,N> &p) const {
            Vec<T,N+1> p_hom(p,1);
            p_hom = inv * p_hom; // this is returning a dynamic matrix for some reason.
            return p_hom.template resized<N>();
        }
        
        /**
         * Inverse transformation of a direction vector; ignores any translation.
         */
        const Vec<T,N> applyInverseVector(const Vec<T,N> &v) const {
            Vec<T,N> o;
            for (index_t r = 0; r < N; r++) {
                for (index_t c = 0; c < N; c++) {
                    o[r] += inv.get(r,c) * v[c];
                }
            }
            return o;
        }
        
        /**
         * Inverse transformation of a normal.
         */
        const Vec<T,N> applyInverseNormal(const Vec<T,N> &n) const {
            // txpose of inverse of inverse.
            Vec<T,N> o;
            for (index_t r = 0; r < N; r++) {
                for (index_t c = 0; c < N; c++) {
                    o[r] += mat.get(c,r) * n[c];
                }
            }
            return o;
        }
        
        /**
         * Concatenation of transforms. 
         * @return A transformation representing a transform by `this` followed by
         * a transform by `at`.
         */
        const AffineTransform<T,N> apply(const AffineTransform<T,N> &at) const {
            AffineTransform atnew;
            atnew.mat = at.mat * mat;
            atnew.inv = inv * at.inv;
            return atnew;
        }
        
        /**
         * @return The inverse transform of `this`.
         */
        const AffineTransform<T,N> inverse() const {
            AffineTransform atnew;
            atnew.mat = inv;
            atnew.inv = mat;
            return atnew;
        }
        
    }; //end AffineTransform class
    
    /*******************************
     * Matrix Construction         *
     *******************************/

    template <typename T> 
    void rotmat(SimpleMatrix<T,4,4> *into, T x, T y, T z, T theta) {
        T c = cos(theta);
        T s = sin(theta);
        T oneMcos = 1 - c;
        T m[4][4] = {
                {c + oneMcos*x*x,           oneMcos*x*y - s*z,      oneMcos*x*z + s*y, 0},
                {oneMcos*y*x + s*z,     c + oneMcos*y*y,            oneMcos*y*z - s*x, 0},
                {oneMcos*z*x - s*y,         oneMcos*z*y + s*x,  c + oneMcos*z*z,       0},
                {0,                         0,                      0,                 1}
        };
        T *from = m[0];
        std::copy(from, from+16, into->begin());
    }
    
    // rotation about an axis
    template <typename T>
    inline void rotmat(SimpleMatrix<T,4,4> *into, const Vec<T,3> &axis, T theta) {
        rotmat(into, axis.x, axis.y, axis.z, theta);
    }

    // rotation from a quaternion
    template <typename T> 
    void rotmat(SimpleMatrix<T,4,4> *into, const Quat<T> &q) {
        T x = q.x; //for convenience
        T y = q.y;
        T z = q.z;
        T w = q.w;
        T x2 = x*x; // for speed
        T y2 = y*y; // ...
        T z2 = z*z;
        T w2 = w*w;
        T twox = 2*x;
        T twoy = 2*y;
        T twoz = 2*z;
        T xy = twox * y;
        T xz = twox * z;
        T xw = twox * w;
        T yz = twoy * z;
        T yw = twoy * w;
        T zw = twoz * w;
        T m[4][4] = {
                {w2 + x2 - y2 - z2, xy - zw, xz + yw, 0},
                {xy + zw, w2 - x2 + y2 - z2, yz - xw, 0},
                {xz - yw, yz + xw, w2 - x2 - y2 + z2, 0},
                {0, 0, 0, 1}
        };
        T *from = m[0];
        std::copy(from, from+16, into->begin());
    }

    // rotation about an arbitrary point
    template <typename T> 
    void rotmat(SimpleMatrix<T,4,4> *into, 
                T x,  T y,  T z, 
                T px, T py, T pz, 
                T theta) {
        T c = cos(theta);
        T s = sin(theta);
        T oneMcos = 1 - c;
        T x2 = x*x;
        T y2 = y*y;
        T z2 = z*z;
        T m[4][4] = {
                {x2 + c*(y2 + z2),   x*y*oneMcos - z*s,  x*z*oneMcos + y*s,  (px*(y2 + z2) - x*(py*y + pz*z))*oneMcos + (py*z - pz*y)*s},
                {x*y*oneMcos + z*s,  y2 + c*(x2 + z2),   y*z*oneMcos - x*s,  (py*(x2 + z2) - y*(px*x + pz*z))*oneMcos + (pz*x - px*z)*s},
                {x*z*oneMcos - y*s,  y*z*oneMcos + x*s,  z2 + c*(x2 + y2),   (pz*(x2 + y2) - z*(px*x + py*y))*oneMcos + (px*y - py*x)*s},
                {0,                  0,                  0,                  1}
        };
        T *from = m[0];
        std::copy(from, from+16, into->begin());
    }
    
    // rotation about an arbitrary point
    template <typename T>
    inline void rotmat(SimpleMatrix<T,4,4> *into,
                       Vec<T,3> axis,
                       Vec<T,3> ctr,
                       T theta) {
        rotmat(into, axis.x, axis.y, axis.z, ctr.x, ctr.y, ctr.z, theta);
    }
    
    // rotation about an arbitrary point (with inverse)
    template <typename T>
    void rotmat(SimpleMatrix<T,4,4> *out_mat, 
                SimpleMatrix<T,4,4> *out_inv, 
                T x,  T y,  T z, 
                T px, T py, T pz, 
                T theta) {
        T c = cos(theta);
        T s = sin(theta);
        T oneMcos = 1 - c;
        T x2 = x*x;
        T y2 = y*y;
        T z2 = z*z;
        T fac[4][4] = {
                {x2 + c*(y2 + z2),   x*y*oneMcos,  x*z*oneMcos,  (px*(y2 + z2) - x*(py*y + pz*z))*oneMcos},
                {x*y*oneMcos,  y2 + c*(x2 + z2),   y*z*oneMcos,  (py*(x2 + z2) - y*(px*x + pz*z))*oneMcos},
                {x*z*oneMcos,  y*z*oneMcos,  z2 + c*(x2 + y2),   (pz*(x2 + y2) - z*(px*x + py*y))*oneMcos},
                {0,                  0,                  0,                  1}
        };
        T sterms[3][4] = {
                {0, -z*s, y*s, (py*z - pz*y)*s},
                {z*s, 0, -x*s, (pz*x - px*z)*s},
                {-y*s, x*s, 0, (px*y - py*x)*s}
        };
        
        // pre-fill output matrices.
        T *from  = fac[0];
        T *sfrom = sterms[0];
        T *mat = out_mat->begin();
        T *inv = out_inv->begin();
        std::copy(from, from+16, mat);
        std::copy(from, from+16, inv);
        
        // add sine terms in.
        for (index_t i = 0; i < 12; i++, mat++, inv++, sfrom++) {
            *mat += *sfrom;
            *inv -= *sfrom; // sin(-x) == -sin(x); thus all the s-terms go negative
        }
    }
    
    // rotation about an arbitrary point (with inverse)
    template <typename T>
    inline void rotmat(SimpleMatrix<T,4,4> *out_mat,
                       SimpleMatrix<T,4,4> *out_inv,
                       const Vec<T,3> &axis,
                       const Vec<T,3> &ctr,
                       T theta) {
        rotmat(out_mat, out_inv, axis.x, axis.y, axis.z, ctr.x, ctr.y, ctr.z, theta);
    }
    
    /*******************************
     * Creation Functions          *
     *******************************/
    
    /** @addtogroup linalg
     *  @{
     */

    /**
     * Rotation about an axis.
     * @param axis Axis of rotation.
     * @param radians Angle of rotation.
     * @return A transformation representing a rotation about `axis` by angle `radians`.
     */
    template <typename T> 
    AffineTransform<T,3> rotation(Vec<T,3> axis, T radians) {
        AffineTransform<T,3> atnew;
        axis = axis.unit();
        rotmat(&atnew.mat, axis.x, axis.y, axis.z, radians);
        transpose(&atnew.inv, atnew.mat); //dst,src. Transpose of a rotation matrix is its inverse
        return atnew;
    }
    
    /**
     * Rotation about point. 
     * 
     * This transformation will not be a pure rotation; it will include a translation
     * component.
     * 
     * @param axis Axis of rotation.
     * @param center Center of rotation.
     * @param radians Angle of rotation.
     * @return A transformation representing a rotation around the point `center` 
     * by angle `radians` and axis `axis`.
     */
    template <typename T> 
    AffineTransform<T,3> rotation(Vec<T,3> axis, const Vec<T,3> &center, T radians) {
        AffineTransform<T,3> atnew;
        axis = axis.unit();
        
        // this is five times faster than translate * rotate * translate
        // and barely faster (3-6%) than calling rotmat() twice with -radians.
        rotmat(&atnew.mat, &atnew.inv, axis.x, axis.y, axis.z, center.x, center.y, center.z, radians);
        return atnew;
    }
    
    /**
     * Rotation from a quaternion.
     * @param q Rotation quaternion.
     * @return A rotation transformation.
     */
    template <typename T> 
    AffineTransform<T,3> rotation(Quat<T> q) {
        AffineTransform<T,3> atnew;
        q = q.unit();
        rotmat(&atnew.mat, q);
        transpose(&atnew.inv, atnew.mat);
        return atnew;
    }
    
    /**
     * 2D rotation about the origin by angle `radians`.
     * @param radians Angle of rotation in the counterclockwise direction
     * @return A 2D rotation transformation.
     */
    template <typename T> 
    AffineTransform<T,2> rotation(T radians) {
        AffineTransform<T,2> atnew;
        T s = sin(radians);
        T c = cos(radians);
        
        atnew.mat.set(0,0, c);
        atnew.mat.set(0,1,-s);
        atnew.mat.set(1,0, s);
        atnew.mat.set(1,1, c);
        
        atnew.inv.set(0,0, c);
        atnew.inv.set(0,1, s);
        atnew.inv.set(1,0,-s);
        atnew.inv.set(1,1, c);
        
        return atnew;
    }
    
    /**
     * Translation transform
     */
    template <typename T, index_t N> 
    AffineTransform<T,N> translation(const Vec<T,N> &tx) {
        AffineTransform<T,N> atnew;
        for (index_t i = 0; i < N; i++) {
            atnew.mat[i][3] =  tx[i];
            atnew.inv[i][3] = -tx[i];
        }
        return atnew;
    }
    
    /**
     * Scale transform. 
     * @param sx Vector whose elements describe a scaling along each axis.
     * @return A transform representing a non-uniform scale along each axis.
     */
    template <typename T, index_t N> 
    AffineTransform<T,N> scale(const Vec<T,N> &sx) {
        AffineTransform<T,N> atnew;
        for (index_t i = 0; i < N; i++) {
            atnew.mat[i][i] = sx[i];
            atnew.inv[i][i] = 1 / sx[i];
        }
        return atnew;
    }
    
    /**
     * Arbitrary transformation
     * @param mat `N x N` matrix representing an arbitrary transformation.
     * @return A transformation by `mat`
     */
    template <typename T, index_t N> 
    AffineTransform<T,N> transformation(const SimpleMatrix<T,N,N> &mat) {
        SimpleMatrix<T,N,N> m_inv;
        AffineTransform<T,N> atnew;
        
        if (N == DYNAMIC_DIM and mat.rows() != mat.cols()) {
            throw NonsquareMatrixException(mat.rows(), mat.cols());
        }
        
        inv(&m_inv, mat);
        MatrixRegion region(MatrixCoord::zeros, MatrixCoord(mat.rows()));
        std::copy(  mat.begin(),   mat.end(), atnew.mat.region(region).first);
        std::copy(m_inv.begin(), m_inv.end(), atnew.inv.region(region).first);
        return atnew;
    }
    
    //2D & 3D convenience:
    
    /**
     * Per-axis 3D scale
     */
    template <typename T> 
    AffineTransform<T,3> scale(T sx, T sy, T sz) {
        return scale(Vec<T,3>(sx,sy,sz));
    }
    
    /**
     * Per-axis 2D scale
     */
    template <typename T> 
    AffineTransform<T,2> scale(T sx, T sy) {
            return scale(Vec<T,2>(sx,sy));
        }
    
    /**
     * 3D translation
     */
    template <typename T> 
    AffineTransform<T,3> translation(T tx, T ty, T tz) {
        return translation(Vec<T,3>(tx,ty,tz));
    }
    
    /**
     * 2D translation
     */
    template <typename T> 
    AffineTransform<T,2> translation(T tx, T ty) {
        return translation(Vec<T,2>(tx,ty));
    }
    
    /// @} // igroup linalg

} //end namespace geom

#endif /* AFFINETRANSFORM_H_ */
