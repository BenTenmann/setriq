//
// Created by Benjamin Tenmann on 18/11/2021.
//

#ifndef METRICS_PAIRWISEDISTANCECOMPUTER_H
#define METRICS_PAIRWISEDISTANCECOMPUTER_H

#include <iostream>
#include <string>
#include <vector>
#include "metrics/Metric.h"
#include "utils/typeDefs.h"

class PairwiseDistanceComputer {
private:
    Metric *distance_metric_;
public:
    PairwiseDistanceComputer() : distance_metric_{} {};

    explicit PairwiseDistanceComputer(Metric* metric) : distance_metric_{metric} {};

    Metric *get_metric() { return this->distance_metric_; };
    double_vector_t compute_distance(const string_vector_t&input_strings) const;
};


#endif //METRICS_PAIRWISEDISTANCECOMPUTER_H
