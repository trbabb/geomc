/*
 * RandomTools.h
 *
 *  Created on: Apr 10, 2011
 *      Author: tbabb
 */

//TODO: test whether using sines/cosines vs. test and discard is faster for 
//      unit ball/sphere random vec generation
//TODO: namespace?

#ifndef RANDOMTOOLS_H_
#define RANDOMTOOLS_H_

#include <vector>

#include <geomc/linalg/Vec.h>
#include <geomc/shape/Rect.h>
#include <geomc/random/Random.h>

namespace geom {

Random* getRandom();

template <typename T> void permute(T objs[], size_t len, Random& rng=*getRandom()){
    for (size_t i = 0; i < len; i++){
        size_t swapIdx = rng.rand(len - i) + i;
        T tmp = objs[swapIdx];
        objs[swapIdx] = objs[i];
        objs[i] = tmp;
    }
}

template <typename T> void permute(std::vector<T> &objs, Random& rng=*getRandom()){
    for (unsigned long i = 0; i < objs.size(); i++){
        unsigned long swapIdx = rng.rand(objs.size() - i) + i;
        T tmp = objs[swapIdx];
        objs[swapIdx] = objs[i];
        objs[i] = tmp;
    }
}

template <typename T>
class RandomVectors {
public:
    Random *rng;
    
    /**********************************
     * Structors                      *
     **********************************/
    
    RandomVectors():rng(getRandom()){
        //do nothing else
    }
    
    RandomVectors(Random *rng):rng(rng){
        //do nothing else
    }
    
    /**********************************
     * Random N-D vectors             *
     **********************************/
    
    // get a random unit vector
    template <index_t N> Vec<T,N> unit(){
        Vec<T,N> v;
        do {
            for (index_t i = 0; i < N; i++){
                // pick a point inside the unit cube
                v[i] = this->rng->template rand<T>(-1,1);
            }
            // discard points outside the unit ball
        } while (v.mag2() > 1);
        return v.unit(); // return point on ball surface
    }
    
    template <index_t N> inline Vec<T,N> unit(T radius){
        return unit<N>() * radius;
    }
    
    template <index_t N> Vec<T,N> solidball(){
        Vec<T,N> v;
        do {
            for (index_t i = 0; i < N; i++){
                v[i] = this->rng->template rand<T>(-1,1);
            }
        } while (v.mag2() > 1);
        return v;
    }
    
    template <index_t N> inline Vec<T,N> solidball(T radius){
        T r = pow(this->rng->template rand<T>(radius), 1.0/N);
        return solidball<N>() * r;
    }
    
    template <index_t N> Vec<T,N> shell(T minradius, T maxradius){
        T radius = pow(this->rng->template rand<T>(minradius, maxradius), 1.0/N);
        return this->unit<N>() * radius;
    }
    
    template <index_t N> Vec<T,N> box(){
        Vec<T,N> v;
        for (index_t i = 0; i < N; i++){
            v[i] = this->rng->template rand<T>(-0.5, 0.5);
        }
        return v;
    }
    
    template <index_t N> Vec<T,N> box(Vec<T,N> lo, Vec<T,N> hi){
        Vec<T,N> v;
        for (index_t i = 0; i < N; i++){  
            v[i] = this->rng->template rand<T>(lo[i], hi[i]);
        }
        return v;
    }
    
    template <index_t N> Vec<T,N> box(const Rect<T,N> &box){
        Vec<T,N> v;
        for (index_t i = 0; i < N; i++){
            v[i] = this->rng->template rand<T>(box.mins[i], box.maxs[i]);
        }
        return v;
    }
    
    /**********************************
     * Random 3D vectors              *
     **********************************/
    
    Vec<T,3> cap(Vec<T,3> dir, T angle, T radius){
        Vec<T,3> v = Vec<T,3>(radius,
                              std::asin(this->rng->template rand<T>(2*angle/M_PI) - 1) + M_PI_2, //phi
                              this->rng->template rand<T>(2*M_PI)).fromSpherical();              //theta
        Vec<T,3> axis = Vec<T,3>(0,0,1).cross(dir);
        if (axis.isZero()) return v;
        return v.rotate(axis, std::acos(dir.z/dir.mag()));
    }
    
    inline Vec<T,3> solidcap(Vec<T,3> dir, T angle, T minR, T maxR){
        T r = this->rng->template rand<T>(minR, maxR);
        r = cbrt(r); //TODO: cbrt before or after?
        return this->cap(dir, angle, r);
    }
    
    // point on a disk of radius <r> with normal <axis>
    Vec<T,3> oriented_disk(Vec<T,3> normal, T r){
        //choose a point on 2D disk first, then
        //rotate the disk to be aligned with <axis>
        Vec<T,3> raxis = Vec<T,3>(0,0,1).cross(normal);
        T angle = this->rng->template rand<T>(2*M_PI);
        T radius = std::sqrt(this->rng->template rand<T>())*r;
        Vec<T,3> disk = this->disk(r).template resized<3>();
        if (raxis.isZero()) return disk;
        return disk.rotate(raxis, std::acos(normal.z/normal.mag()));
    }

    inline Vec<T,3> box(T minx, T miny, T minz,
                        T maxx, T maxy, T maxz){
        return Vec<T,3>(
                this->rng->template rand<T>(minx, maxx),
                this->rng->template rand<T>(miny, maxy),
                this->rng->template rand<T>(minz, maxz)
        );
    }
    
    /**********************************
     * Random 2D vectors              *
     **********************************/
    
    inline Vec<T,2> disk(){
        return solidball<2>();
    }
    
    inline Vec<T,2> disk(T radius){
        return solidball<2>() * radius;
    }
    
    inline Vec<T,2> arc(T min_radians, T max_radians){
        T angle = this->rng->template rand<T>(min_radians, max_radians);
        return Vec<T,2>(-cos(angle), sin(angle));
    }
    
    Vec<T,2> solidarc(T min_radians, T max_radians,
                      T min_radius,  T max_radius){
        T angle = this->rng->template rand<T>(min_radians, max_radians);
        T len   = std::sqrt(this->rng->template rand<T>(min_radius, max_radius));
        return Vec<T,2>(len,angle).fromPolar();
    }
    
    inline Vec<T,2> box(T minx, T miny, 
                        T maxx, T maxy){
        return Vec<T,2>(this->rng->template rand<T>(minx, maxx),
                        this->rng->template rand<T>(miny, maxy));
    }
};

}; //end namespace geom

#endif /* RANDOMTOOLS_H_ */
