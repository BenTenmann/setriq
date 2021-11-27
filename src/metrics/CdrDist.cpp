//
// Created by Benjamin Tenmann on 20/11/2021.
//

#include <cmath>
#include <utility>
#include "metrics/CdrDist.h"

metric::CdrDist::CdrDist(const doubleMatrix& matrix, const stringIndexMap& index, double gapPenalty, size_t cacheSize) {
    SubstitutionMatrix substitutionMatrix(matrix, index);
    SmithWaterman algo(substitutionMatrix, gapPenalty, cacheSize);

    this->algorithm = algo;
}

double metric::CdrDist::forward(const std::string &a, const std::string &b) {
    double abScore = this->algorithm(a, b);
    double aaScore = this->algorithm.identityScore(a);
    double bbScore = this->algorithm.identityScore(b);

    return 1 - std::sqrt((abScore * abScore) / (aaScore * bbScore));
}
