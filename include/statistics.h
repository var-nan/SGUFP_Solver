//
// Created by nandgate on 4/1/2025.
//

#ifndef STATISTICS_H
#define STATISTICS_H
#include <algorithm>

static double stats_mean(const double *data, size_t size){
    double sum = 0;
    for (int i = 0; i < size; i++) {
        sum += data[i];
    }
    return sum/size;
}

static double stats_min(const double *data, size_t size) {
    double min_val = data[0];
    for (int i = 1; i < size; i++)
        min_val = std::min(min_val, data[i]);
    return min_val;
}

static double stats_max(const double *data, size_t size) {
    double max_val = data[0];
    for (int i = 1; i < size; i++)
        max_val = std::max(max_val, data[i]);
    return max_val;
}

static double stats_absdev(const double *data, size_t size) {
    return 0;
}

#endif //STATISTICS_H
