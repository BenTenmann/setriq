//
// Created by Benjamin Tenmann on 20/11/2021.
//

#include "alignment/SubstitutionMatrix.h"

SubstitutionMatrix::SubstitutionMatrix(const double_matrix_t& matrix, const token_index_map_t& index) {
    /**
     * Initialize a SubstitutionMatrix object.
     *
     * @param matrix: the substitution scoring matrix
     * @param index: the token index map
     */
    this->scoring_matrix_ = matrix;
    this->token_map_ = index;
}

double SubstitutionMatrix::forward(const char &from, const char &to) {
    /**
     * Retrieve the substitution score for two input characters.
     *
     * @param from: the first character in the substitution
     * @param to: the second character in the substitution
     * @return the substitution score between two characters
     */
    size_t fromIdx, toIdx;
    fromIdx = this->token_map_[from];
    toIdx = this->token_map_[to];

    return this->scoring_matrix_[fromIdx][toIdx];
}
