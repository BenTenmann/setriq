//
// Created by Benjamin Tenmann on 18/11/2021.
//

#include "PairwiseDistanceComputer.h"

PairwiseDistanceComputer::PairwiseDistanceComputer(Metric* metric) {
    this->distanceMetric = metric;
}

doubleVector PairwiseDistanceComputer::computeDistance(const stringVector& inputStrings) const {
    /**
     * Method for computing a given metric on a set of sequences, pairwise. In case OpenMP is available at compile-time,
     * the method will be multi-threaded.
     *
     * @param inputStrings: a vector of input strings
     * @return a flat (N * (N - 1)) / 2 vector of doubles
     */
    size_t N = inputStrings.size();
    doubleVector distanceMatrix(N * (N - 1) / 2);

    #pragma omp parallel for default(none) shared(inputStrings, distanceMatrix, N)
    for (size_t i = 0; i < N - 1; i++) {
        for (size_t j = i + 1; j < N; j++) {
            size_t index = (N * (N - 1)) / 2 - (N - i) * ((N - i) - 1) / 2 + j - i - 1;
            distanceMatrix[index] = (*this->distanceMetric)(inputStrings[i], inputStrings[j]);
        }
    }

    return distanceMatrix;
}
