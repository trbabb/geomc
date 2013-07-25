/*
 * RandomBattery.h
 *
 *  Created on: Nov 11, 2012
 *      Author: tbabb
 */

#ifndef RANDOMBATTERY_H_
#define RANDOMBATTERY_H_

#include <geomc/random/Random.h>

namespace geom {

struct result_t {
    double stdev;
    double mean;
};

class RandomBattery {
public:
    RandomBattery(Random *rng);
    virtual ~RandomBattery();
    
    result_t distribute(int trials);
    result_t bitHistogram(int trials);
    result_t byteHistogram(int trials);
    double speed(int trials);
    
    void runAll(int trials);
    
    void plot(double *bins, int n_bins, double max_value);
    
protected:
    Random *rng;
};

} /* namespace geom */
#endif /* RANDOMBATTERY_H_ */
