//
// Created by Benjamin Tenmann on 21/11/2021.
//


#include <utility>

#include "alignment/SmithWaterman.h"

SmithWaterman::SmithWaterman(SubstitutionMatrix matrix, const double& gap_opening_penalty, const double& gap_extension_penalty)
    : gap_opening_penalty_{gap_opening_penalty}, gap_extension_penalty_{gap_extension_penalty} {
    /**
     * Initialize a SmithWaterman object.
     *
     * @param matrix: a SubstitutionMatrix objet defining the substitution scores in the alignment
     * @param gap_penalty: the penalty for a gap in the alignment
     */
    this->substitution_matrix_ = std::move(matrix);
}

double SmithWaterman::fill_scoring_matrix_(const std::string &a,
                                           const std::string &b) const {
    /**
     * Fill the alignment scoring matrix and return the maximal alignment score between two input strings.
     *
     * @param a: an input string to be aligned
     * @param b: an input string to be aligned
     * @return the maximal alignment score between two sequences
     */
    const auto& n = a.size();
    const auto& m = b.size();

    constexpr auto k_axis = 0ul;
    constexpr auto l_axis = 1ul;

    auto* ptr_a = &a.front();
    auto* ptr_b = &b.front();

    auto&& max_score = 0.;
    auto&& scoring_matrix = double_matrix_t (n + 1, double_vector_t (m + 1, 0));
    for (size_t i = 1; i < (n + 1); i++) {
        for (size_t j = 1; j < (m + 1); j++) {
            const auto& _i = i - 1;
            const auto& _j = j - 1;

            const auto& alignment_score = scoring_matrix[_i][_j] + this->substitution_matrix_(*(ptr_a + _i), *(ptr_b + _j));
            const auto& k_gap_score = this->calculate_gap_penalty_(scoring_matrix, i, j, k_axis);
            const auto& l_gap_score = this->calculate_gap_penalty_(scoring_matrix, j, i, l_axis);

            const auto& current_score = std::max(alignment_score, std::max(k_gap_score, l_gap_score));
            if (current_score > max_score)
                max_score = current_score;

            scoring_matrix[i][j] = current_score;
        }
    }
    return max_score;
}

double SmithWaterman::calculate_gap_penalty_(const double_matrix_t &scoring_matrix,
                                             const size_t &max_gap_length,
                                             const size_t &idx,
                                             const size_t &axis) const {
    /**
     * Calculate the gap score along a give axis. It has a lower bound of 0 and is linear (non-affine).
     *
     * @param scoring_matrix: the alignment scoring matrix
     * @param max_gap_length: the maximum possible gap length
     * @param idx: the current idx position (row or column)
     * @param axis: the axis over which to iterate
     * @return the maximum gap-score
     */
    auto&& max_score = 0.;
    for (size_t i = 1; i < (max_gap_length + 1); i++) {
        const auto& k = max_gap_length - i;
        const auto& elem = axis ? scoring_matrix[idx][k] : scoring_matrix[k][idx];
        const auto& score = elem - this->gap_opening_penalty_ - ((double)(i - 1) * this->gap_extension_penalty_);

        if (score > max_score)
            max_score = score;
    }
    return max_score;
}

double SmithWaterman::compute_best_alignment_score_(const std::string &a, const std::string &b) const {
    /**
     * Compute the maximal alignment score between two strings
     *
     * @param a: an input string to be aligned
     * @param b: an input string to be aligned
     * @return the maximal alignment score
     */
    auto&& score = this->fill_scoring_matrix_(a, b);

    return score;
}

double SmithWaterman::forward(const std::string &a, const std::string &b) const {
    /**
     * Wrapper for computing the maximal alignment score
     *
     * @param a: an input string to be aligned
     * @param b: an input string to be aligned
     * @return the maximal alignment score
     */
    const auto& best_score = this->compute_best_alignment_score_(a, b);
    return best_score;
}

double SmithWaterman::identity_score(const std::string &input_string) const {
    /**
     * Compute the maximal alignment score for an input string with itself. It is a trivial case of SW, as the maximal
     * score will always be at the final diagonal position of the scoring matrix. Thus it reduces down to a cumulative
     * sum of the substitution scores of the string characters with themselves.
     *
     * @param input_string: an input string to be aligned with itself
     * @return the maximal self-alignment score for an input string
     */
    const auto* end = &input_string.back() + 1;

    auto&& score = 0.;
    for (auto* ptr = &input_string.front(); ptr != end; ptr++) {
        score += this->substitution_matrix_(*ptr, *ptr);
    }
    return score;
}
