/*
 * RandomBattery.cpp
 *
 *  Created on: Nov 11, 2012
 *      Author: tbabb
 */

#include "RandomBattery.h"
#include <math.h>
#include <time.h>
#include <iostream>

namespace geom {

void stdev(result_t &result, double *nums, int n_nums){
    double mean = 0;
    for (int i = 0; i < n_nums; i++){
        mean += nums[i];
    }
    mean /= n_nums;
    
    double stdev_sum = 0;
    for (int i = 0; i < n_nums; i++){
        double delta = nums[i] - mean; 
        stdev_sum += delta * delta;
    }
    
    result.stdev = sqrt(stdev_sum / n_nums);
    result.mean = mean;
    
}

RandomBattery::RandomBattery(Random *rng):rng(rng) {}

RandomBattery::~RandomBattery() {}

result_t RandomBattery::distribute(int trials){
    const int n_bins = 32;
    double bins[n_bins];
    
    for (int i = 0; i < n_bins; i++){
        bins[i] = 0;
    }
    
    for (int i = 0; i < trials; i++){
        int n = rng->rand<int>(n_bins);
        bins[n]++;
    }
    
    plot(bins, n_bins, (float)trials/n_bins);
    result_t result;
    stdev(result, bins, n_bins);
    return result;
}

result_t RandomBattery::bitHistogram(int trials){
    double places[32];
    for (int i = 0; i < 32; i++){
        places[i] = 0;
    }
    for (int i = 0; i < trials; i++){
        uint32_t x = rng->rand32();
        for (int b = 0; b < 32; b++){
            places[b] += x & 1;
            x = x >> 1;
        }
    }
    
    plot(places, 32, trials/2);
    result_t result;
    stdev(result, places, 32);
    return result;
}

result_t RandomBattery::byteHistogram(int trials){
    double bytes[256];
    for (int i = 0; i < 256; i++){
        bytes[i] = 0;
    }
    
    for (int i = 0; i < trials; i++){
        uint32_t x = rng->rand32();
        bytes[x & 0x000000ff      ]++;
        bytes[x & 0x0000ff00 >>  8]++;
        bytes[x & 0x00ff0000 >> 16]++;
        bytes[x & 0xff000000 >> 24]++;
    }
    
    plot(bytes, 256, 4*trials/256);
    result_t result;
    stdev(result, bytes, 256);
    return result;
}

double RandomBattery::speed(int trials){
    uint32_t intermediate = 0;
    clock_t begin = clock();
    for (int i = 0; i < trials; i++){
        intermediate = intermediate ^ rng->rand32();
    }
    clock_t end = clock();
    return (end-begin) / (double)CLOCKS_PER_SEC;
}

void RandomBattery::runAll(int trials){
    result_t rslt;
    
    std::cout << "distribution over N buckets:" << std::endl;
    rslt = distribute(trials);
    std::cout << "stdev: " << rslt.stdev << " mean: " << rslt.mean << " spread: " << (rslt.stdev / rslt.mean) << std::endl;

    std::cout << "distribution of bits:" << std::endl;
    rslt = bitHistogram(trials);
    std::cout << "stdev: " << rslt.stdev << " mean: " << rslt.mean << " spread: " << (rslt.stdev / rslt.mean) << std::endl;

    std::cout << "distribution of bytes:" << std::endl;
    rslt = byteHistogram(trials);
    std::cout << "stdev: " << rslt.stdev << " mean: " << rslt.mean << " spread: " << (rslt.stdev / rslt.mean) << std::endl;
    
    double t = speed(trials);
    std::cout << "time: " << t << " (" << (trials/t) << " words/sec)" << std::endl;
}

void RandomBattery::plot(double *bins, int n_bins, double max_value){
    int char_width = 80;
    for (int bin = 0; bin < n_bins; bin++){
        int stars = round(char_width * (bins[bin] / max_value));
        for (int s = 0; s < stars; s++){
            std::cout << "*";
        }
        std::cout << std::endl;
    }
}



} /* namespace geom */
