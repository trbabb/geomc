/*
 * MTRand.cpp
 * Mersenne Twister Random Number Generator
 *
 *  Created on: Feb 22, 2009
 *      Author: Tim Babb
 *
 * Based on public domain code by michael brundage. snipped from:
 * http://www.qbrundage.com/michaelb/pubs/essays/random_number_generation.html
 * 
 * NOT CRYPTOGRAPHICALLY SECURE! Don't use for cryptography/gambling
 */

#include <cstdlib>
#include <ctime>
#include <limits>
#include "random/MTRand.h"
#include "random/LCRand.h"

using namespace geom;

int MTRand::uniq = 0;

#define MT_IA           397
#define MT_IB           (MT_LEN - MT_IA)
#define UPPER_MASK      0x80000000
#define LOWER_MASK      0x7FFFFFFF
#define MATRIX_A        0x9908B0DF
#define TWIST(b,i,j)    ((b)[i] & UPPER_MASK) | ((b)[j] & LOWER_MASK)
#define MAGIC(s)        (((s)&1)*MATRIX_A)

MTRand::MTRand(){
    this->rseed(time(0) ^ (clock() + uniq++));
}

MTRand::MTRand(uint64_t seed){
    this->rseed(seed);
}

MTRand::~MTRand() {
    //do nothing
}

uint32_t MTRand::rand32() {
    uint32_t *b = mt_buffer;
    int idx = mt_index;
    uint32_t s;
    int i;

    if (idx == MT_LEN*sizeof(uint32_t)) {
        idx = 0;
        i = 0;
        for (; i < MT_IB; i++) {
            s = TWIST(b, i, i+1);
            b[i] = b[i + MT_IA] ^ (s >> 1) ^ MAGIC(s);
        }
        for (; i < MT_LEN-1; i++) {
            s = TWIST(b, i, i+1);
            b[i] = b[i - MT_IB] ^ (s >> 1) ^ MAGIC(s);
        }

        s = TWIST(b, MT_LEN-1, 0);
        b[MT_LEN-1] = b[MT_IA-1] ^ (s >> 1) ^ MAGIC(s);
    }
    mt_index = idx + sizeof(uint32_t);
    return *(uint32_t*)((unsigned char*)b + idx);
}

void MTRand::rseed(uint64_t seed){
    LCRand seeder(seed);
    for (int i = 0; i < MT_LEN; i++){
        mt_buffer[i] = seeder.rand32();
    }
    mt_index = 0;
}
