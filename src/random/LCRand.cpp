/*
 * LCRand.cpp
 *
 *  Created on: Nov 11, 2012
 *      Author: tbabb
 */

#include <ctime>
#include "random/LCRand.h"

namespace geom {

int LCRand::uniq = 0;

LCRand::LCRand():state(std::time(0) ^ (std::clock() + uniq++)) {}

LCRand::LCRand(uint64_t seed):state(seed) {}

LCRand::~LCRand() {/* do nothing */}

uint32_t LCRand::rand32(){
    return rand64() >> 32; // return high quality/period bits
}

void LCRand::rseed(uint64_t seed){
    state = seed;
}

} /* namespace geom */
