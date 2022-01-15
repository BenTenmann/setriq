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

double SmithWaterman::fill_scoring_matrix_(double_matrix_t &scoring_matrix,
                                           const std::string &a,
                                           const std::string &b) {
    /**
     * Fill the alignment scoring matrix and return the maximal alignment score between two input strings.
     *
     * @param scoring_matrix: the unfilled scoring matrix
     * @param a: an input string to be aligned
     * @param b: an input string to be aligned
     * @return the maximal alignment score between two sequences
     */
    const auto& N = scoring_matrix.size();
    const auto& M = scoring_matrix[0].size();

    double alignment_score, k_gap_score, l_gap_score;
    double current_score, max_score {0};

    const size_t& k_axis = 0;
    const size_t& l_axis = 1;

    for (size_t i = 1; i < N; i++) {
        for (size_t j = 1; j < M; j++) {
            alignment_score = scoring_matrix[i - 1][j - 1] + this->substitution_matrix_(a[i - 1], b[j - 1]);
            k_gap_score = this->calculate_gap_penalty_(scoring_matrix, i, j, k_axis);
            l_gap_score = this->calculate_gap_penalty_(scoring_matrix, j, i, l_axis);

            current_score = std::max(alignment_score, std::max(k_gap_score, l_gap_score));
            if (current_score > max_score)
                max_score = current_score;

            scoring_matrix[i][j] = current_score;
        }
    }
    return max_score;
}

double SmithWaterman::calculate_gap_penalty_(const double_matrix_t &scoring_matrix,
                                             const size_t &max_gap_length,
                                             const size_t &index,
                                             const size_t &axis) const {
    /**
     * Calculate the gap score along a give axis. It has a lower bound of 0 and is linear (non-affine).
     *
     * @param scoring_matrix: the alignment scoring matrix
     * @param max_gap_length: the maximum possible gap length
     * @param index: the current index position (row or column)
     * @param axis: the axis over which to iterate
     * @return the maximum gap-score
     */
    double max_score {0};

    double elem, score;
    size_t k;
    for (size_t i = 1; i < (max_gap_length + 1); i++) {
        k = max_gap_length - i;
        elem = axis ? scoring_matrix[index][k] : scoring_matrix[k][index];
        score = elem - this->gap_opening_penalty_ - ((double)(i - 1) * this->gap_extension_penalty_);

        if (score > max_score)
            max_score = score;
    }
    return max_score;
}

double SmithWaterman::compute_best_alignment_score_(const std::string &a, const std::string &b) {
    /**
     * Compute the maximal alignment score between two strings
     *
     * @param a: an input string to be aligned
     * @param b: an input string to be aligned
     * @return the maximal alignment score
     */
    const auto& N = a.size();
    const auto& M = b.size();

    double_matrix_t scoring_matrix (N + 1, double_vector_t (M + 1, 0));
    const double& score {this->fill_scoring_matrix_(scoring_matrix, a, b)};

    return score;
}

double SmithWaterman::forward(const std::string &a, const std::string &b) {
    /**
     * Wrapper for computing the maximal alignment score
     *
     * @param a: an input string to be aligned
     * @param b: an input string to be aligned
     * @return the maximal alignment score
     */
    const auto& best_score {this->compute_best_alignment_score_(a, b)};
    return best_score;
}

double SmithWaterman::identity_score(const std::string &input_string) {
    /**
     * Compute the maximal alignment score for an input string with itself. It is a trivial case of SW, as the maximal
     * score will always be at the final diagonal position of the scoring matrix. Thus it reduces down to a cumulative
     * sum of the substitution scores of the string characters with themselves.
     *
     * @param input_string: an input string to be aligned with itself
     * @return the maximal self-alignment score for an input string
     */
    const auto& N = input_string.size();

    double score {0};
    for (size_t i = 0; i < N; i++) {
        score += this->substitution_matrix_(input_string[i], input_string[i]);
    }
    return score;
}
