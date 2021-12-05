//
// Created by Benjamin Tenmann on 20/11/2021.
//

#ifndef METRICS_SUBSTITUTIONMATRIX_H
#define METRICS_SUBSTITUTIONMATRIX_H

#include <string>
#include "utils/typeDefs.h"

class SubstitutionMatrix {
private:
    doubleMatrix subMatrix;
    stringIndexMap edgeMap;

public:
    SubstitutionMatrix() : subMatrix(), edgeMap() {};
    SubstitutionMatrix(const doubleMatrix&, const stringIndexMap&);

    double forward (const char&, const char&);
    double operator () (const char& a, const char &b) { return this->forward(a, b); };
};


#endif //METRICS_SUBSTITUTIONMATRIX_H
