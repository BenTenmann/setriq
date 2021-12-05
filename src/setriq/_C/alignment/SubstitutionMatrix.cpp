//
// Created by Benjamin Tenmann on 20/11/2021.
//

#include "alignment/SubstitutionMatrix.h"

SubstitutionMatrix::SubstitutionMatrix(const doubleMatrix& matrix, const stringIndexMap& index) {
    this->subMatrix = matrix;
    this->edgeMap = index;
}

double SubstitutionMatrix::forward(const char &from, const char &to) {
    size_t fromIdx, toIdx;
    fromIdx = this->edgeMap[from];
    toIdx = this->edgeMap[to];

    return this->subMatrix[fromIdx][toIdx];
}
