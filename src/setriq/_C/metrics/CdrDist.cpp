//
// Created by Benjamin Tenmann on 20/11/2021.
//

#include <cmath>
#include <utility>
#include "metrics/CdrDist.h"

metric::CdrDist::CdrDist(const doubleMatrix& matrix, const stringIndexMap& index, double gapPenalty) {
    /**
     * Initialize a CdrDist object. Uses the SmithWaterman algorithm for the sequence alignment.
     *
     * @param matrix: the substitution scoring matrix
     * @param index: the token-index map
     * @param gapPenalty: the penalty for a gap in the alignment
     */
    SubstitutionMatrix substitutionMatrix(matrix, index);
    SmithWaterman algo(substitutionMatrix, gapPenalty);

    this->algorithm = algo;
}

double metric::CdrDist::forward(const std::string &a, const std::string &b) {
    /**
     * Compute the CdrDist metric between two input strings.
     *
     * @param a: on input string to be compared
     * @param b: second input string to be compared
     * @return the CdrDist metric between the input strings
     */

    // calculate the numerator of CdrDist. This is the most expensive operation of the metric, as the sequences get
    // aligned
    double abScore = this->algorithm(a, b);

    // when a == b in SW, the best score collapses down to a simple arithmetic sum over all the amino acid positions
    // this gives a significant speed bump
    double aaScore = this->algorithm.identityScore(a);
    double bbScore = this->algorithm.identityScore(b);

    return 1 - std::sqrt((abScore * abScore) / (aaScore * bbScore));
}
