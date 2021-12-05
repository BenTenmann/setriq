//
// Created by Benjamin Tenmann on 19/11/2021.
//

#include <stdexcept>
#include "metrics/TcrDist.h"

metric::TcrDist::TcrDist(const doubleMatrix& subMat,
                         const stringIndexMap& index,
                         double gapPen,
                         char gapSym,
                         double weight) {
    /**
     * Initialize the TcrDist object.
     *
     * @param subMat: the substitution scoring matrix
     * @param index: the token index
     * @param gapPen: the metric gap penalty
     * @param gapSym: the gap symbol to be used (e.g. "-")
     * @param weight: the weight of the metric component output
     */

    SubstitutionMatrix sm {subMat, index};
    this->substitutionMatrix = sm;
    this->gapPenalty = gapPen;
    this->gapSymbol = gapSym;
    this->distanceWeight = weight;
}

double metric::TcrDist::forward(const std::string &a, const std::string &b) {
    /**
     * Compute the TcrDist metric (Dash et al) between two provided sequences. Be sure to provide sequences of equal
     * length. Errors are handled in the Python API for speed and interpretability considerations.
     *
     * @param a: a string to be compared
     * @param b: another string to be compared
     * @return TcrDist metric between the two strings
     */

    double distance {0}, substitution;
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i] == b[i]) continue;

        if (a[i] == this->gapSymbol || b[i] == this->gapSymbol) {
            distance += this->gapPenalty;
            continue;
        }

        substitution = this->substitutionMatrix(a[i], b[i]);
        distance += std::min(4., 4. - substitution);
    }
    return distance * this->distanceWeight;
}
