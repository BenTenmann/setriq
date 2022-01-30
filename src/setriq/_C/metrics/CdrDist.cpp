//
// Created by Benjamin Tenmann on 20/11/2021.
//

#include <cmath>
#include <utility>

#include "metrics/CdrDist.h"

metric::CdrDist::CdrDist(const double_matrix_t& matrix,
                         const token_index_map_t& index,
                         const double& gap_opening_penalty,
                         const double& gap_extension_penalty) {
    /**
     * Initialize a CdrDist object. Uses the SmithWaterman algorithm_ for the sequence alignment.
     *
     * @param matrix: the substitution scoring matrix
     * @param index: the token-index map
     * @param gap_penalty: the penalty for a gap in the alignment
     */
    SubstitutionMatrix substitution_matrix {matrix, index};
    this->algorithm_ = SmithWaterman (substitution_matrix, gap_opening_penalty, gap_extension_penalty);
}

double metric::CdrDist::forward(const std::string &a, const std::string &b) const {
    /**
     * Compute the CdrDist metric between two input strings.
     *
     * @param a: on input string to be compared
     * @param b: second input string to be compared
     * @return the CdrDist metric between the input strings
     */

    // calculate the numerator of CdrDist. This is the most expensive operation of the metric, as the sequences get
    // aligned
    const auto& ab_score = this->algorithm_(a, b);

    // when a == b in SW, the best score collapses down to a simple arithmetic sum over all the amino acid positions
    // this gives a significant speed bump
    const auto& aa_score = this->algorithm_.identity_score(a);
    const auto& bb_score = this->algorithm_.identity_score(b);

    return 1 - std::sqrt((ab_score * ab_score) / (aa_score * bb_score));
}
