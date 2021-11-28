//
// Created by Benjamin Tenmann on 19/11/2021.
//

#include <stdexcept>
#include "metrics/TcrDist.h"

metric::TcrDist::TcrDist(const doubleMatrix& subMat, const stringIndexMap& index, double gapPen, char gapSym, double weight) {
    SubstitutionMatrix sm {subMat, index};
    this->substitutionMatrix = sm;
    this->gapPenalty = gapPen;
    this->gapSymbol = gapSym;
    this->distanceWeight = weight;
}

double metric::TcrDist::forward(const std::string &a, const std::string &b) {
    if (a.size() != b.size()) throw std::invalid_argument("input strings must be of equal length.");

    double distance {0}, substitution;
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i] == b[i]) continue;

        if (a[i] == this->gapSymbol || b[i] == this->gapSymbol) {
            distance += this->gapPenalty * this->distanceWeight;
            continue;
        }

        substitution = this->substitutionMatrix(a.substr(i, 1), b.substr(i, 1));
        distance += std::min(4., 4. - substitution) * this->distanceWeight;
    }
    return distance;
}
