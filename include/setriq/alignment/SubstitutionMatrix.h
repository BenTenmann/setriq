//
// Created by Benjamin Tenmann on 20/11/2021.
//

#ifndef METRICS_SUBSTITUTIONMATRIX_H
#define METRICS_SUBSTITUTIONMATRIX_H

#include <string>

#include "utils/type_defs.h"

class SubstitutionMatrix {
private:
    double_matrix_t scoring_matrix_;
    token_index_map_t token_map_;

public:
    SubstitutionMatrix() : scoring_matrix_{}, token_map_{} {};
    SubstitutionMatrix(const double_matrix_t&, const token_index_map_t&);

    double forward (const char&, const char&) const;
    double operator () (const char& a, const char& b) const { return this->forward(a, b); };
};


#endif //METRICS_SUBSTITUTIONMATRIX_H
