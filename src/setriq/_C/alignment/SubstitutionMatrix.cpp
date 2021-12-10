//
// Created by Benjamin Tenmann on 20/11/2021.
//

#include "alignment/SubstitutionMatrix.h"

SubstitutionMatrix::SubstitutionMatrix(const doubleMatrix& matrix, const stringIndexMap& index) {
    /**
     * Initialize a SubstitutionMatrix object.
     *
     * @param matrix: the substitution scoring matrix
     * @param index: the token index map
     */
    this->subMatrix = matrix;
    this->edgeMap = index;
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
    fromIdx = this->edgeMap[from];
    toIdx = this->edgeMap[to];

    return this->subMatrix[fromIdx][toIdx];
}
