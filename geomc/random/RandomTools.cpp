/*
 * RandomTools.cpp
 *
 *  Created on: Apr 10, 2011
 *      Author: tbabb
 */

#include <cmath>
#include <geomc/random/RandomTools.h>
#include <geomc/random/MTRand.h>

using namespace geom;

//lazily init'd from getRandom()
Random *randomtools_rng = NULL;

Random* geom::getRandom(){
    if (randomtools_rng == NULL){
        randomtools_rng = new MTRand();
    }
    return randomtools_rng;
}
