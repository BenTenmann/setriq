//
// Created by Benjamin Tenmann on 18/11/2021.
//

#include "PairwiseDistanceComputer.h"

PairwiseDistanceComputer::PairwiseDistanceComputer(Metric* metric) {
    this->distanceMetric = metric;
}

doubleVector PairwiseDistanceComputer::computeDistance(const stringVector& inputStrings) const {
    size_t N = inputStrings.size();
    doubleVector distanceMatrix(N * (N - 1) / 2);

    double *ptr = &distanceMatrix.front();
    for (size_t i = 0; i < N - 1; i++) {
        for (size_t j = i + 1; j < N; j++) {
            *ptr = (*this->distanceMetric)(inputStrings[i], inputStrings[j]);
            ptr++;
        }
    }

    return distanceMatrix;
}
