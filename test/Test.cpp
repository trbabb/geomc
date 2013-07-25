/*
 * Test.cpp
 *
 *  Created on: Nov 21, 2010
 *      Author: tbabb
 */

#include <iostream>

#include <geomc/function/Raster.h>

using namespace std;
using namespace geom;

typedef Raster<double,float,2,3> Image2d3f;

int main_test(int arc, char** argv){
//int main(int arc, char** argv){
    Vec<size_t,2> d = Vec<size_t,2>(512,256);
    Image2d3f img = Image2d3f(d);
    
    img.set(Vec2i(10,10), Vec3f(1,0,0));
    
    
    for (int i = 0; i < 40; i++){
        double s = 8 + i/10.0;
        cout << s << " " << img.sample<EDGE_CLAMP,INTERP_LINEAR>(Vec2d(10.62314,s)) << endl;
    }
    
    double b[4][4] = {
            {0.10, 0.18, 0.31, 0.20},
            {0.16, 0.23, 0.52, 0.35},
            {0.30, 0.50, 0.57, 0.21},
            {0.12, 0.36, 0.42, 0.14}
    };
    
    cout << endl << "---------" << endl << endl;
    
    for (int i = 0; i < 10; i++){
        double s[2] = {i / 10.0, 0.5}; 
        cout << s[0] << " " << interp_cubic(s, b[0], 2) << std::endl;
    }
    
    return 0;
}
