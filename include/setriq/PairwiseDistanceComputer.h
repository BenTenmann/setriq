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
    Metric *distanceMetric;
public:
    PairwiseDistanceComputer() : distanceMetric() {};

    explicit PairwiseDistanceComputer(Metric*);

    Metric *getMetric() { return this->distanceMetric; };
    doubleVector computeDistance(const stringVector&) const;
};


#endif //METRICS_PAIRWISEDISTANCECOMPUTER_H
