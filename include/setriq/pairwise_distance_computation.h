//
// Created by Benjamin Tenmann on 29/01/2022.
//

#ifndef SETRIQ_PAIRWISE_DISTANCE_COMPUTATION_H
#define SETRIQ_PAIRWISE_DISTANCE_COMPUTATION_H

#include "utils/type_defs.h"

template<typename T>
double_vector_t pairwise_distance_computation(T metric, const string_vector_t& input_strings) {
    const auto& n = input_strings.size();
    auto&& distance_matrix = double_vector_t (n * (n - 1) / 2);

#pragma omp parallel for default(none) firstprivate(n, metric) shared(input_strings, distance_matrix)
    for (size_t i = 0; i < (n - 1); i++) {
        for (size_t j = (i + 1); j < n; j++) {
            const auto& idx = (n * (n - 1)) / 2 - (n - i) * ((n - i) - 1) / 2 + j - i - 1;
            distance_matrix[idx] = metric.forward(input_strings[i], input_strings[j]);
        }
    }
    return distance_matrix;
}

#endif //SETRIQ_PAIRWISE_DISTANCE_COMPUTATION_H
