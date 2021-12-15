//
// Created by Benjamin Tenmann on 19/11/2021.
//

#include <stdexcept>
#include "metrics/TcrDist.h"

metric::TcrDist::TcrDist(const double_matrix_t& scoring_matrix,
                         const token_index_map_t& index,
                         double gap_penalty,
                         char gap_symbol,
                         double weight) : gap_penalty_{gap_penalty}, gap_symbol_{gap_symbol}, distance_weight_{weight} {
    /**
     * Initialize the TcrDist object.
     *
     * @param scoring_matrix: the substitution scoring matrix
     * @param index: the token index
     * @param gap_penalty: the metric gap penalty
     * @param gap_symbol: the gap symbol to be used (e.g. "-")
     * @param weight: the weight of the metric component output
     */
    this->substitution_matrix_ = SubstitutionMatrix (scoring_matrix, index);
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
    const auto& max_distance {4.};

    double distance {0}, substitution;
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i] == b[i]) continue;

        if (a[i] == this->gap_symbol_ || b[i] == this->gap_symbol_) {
            distance += this->gap_penalty_;
            continue;
        }

        substitution = max_distance - this->substitution_matrix_(a[i], b[i]);
        distance += std::min(max_distance, substitution);
    }
    return distance * this->distance_weight_;
}
