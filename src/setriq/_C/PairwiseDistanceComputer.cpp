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

    #pragma omp parallel for default(none) shared(inputStrings, distanceMatrix, N, std::cout)
    for (size_t i = 0; i < N - 1; i++) {
        for (size_t j = i + 1; j < N; j++) {
            size_t index = (N * (N - 1)) / 2 - (N - i) * ((N - i) - 1) / 2 + j - i - 1;
            distanceMatrix[index] = (*this->distanceMetric)(inputStrings[i], inputStrings[j]);
        }
    }

    return distanceMatrix;
}
