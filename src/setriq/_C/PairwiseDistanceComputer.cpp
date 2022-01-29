//
// Created by Benjamin Tenmann on 18/11/2021.
//

#include "PairwiseDistanceComputer.h"

double_vector_t PairwiseDistanceComputer::compute_distance(const string_vector_t& input_strings) const {
    /**
     * Method for computing a given metric on a set of sequences, pairwise. In case OpenMP is available at compile-time,
     * the method will be multi-threaded.
     *
     * @param input_strings: a vector of input strings
     * @return a flat (n * (n - 1)) / 2 vector of doubles
     */
    const auto& n = input_strings.size();
    double_vector_t distance_matrix(n * (n - 1) / 2);

#pragma omp parallel for default(none) firstprivate(n) shared(input_strings, distance_matrix)
    for (size_t i = 0; i < n - 1; i++) {
        for (size_t j = i + 1; j < n; j++) {
            size_t index = (n * (n - 1)) / 2 - (n - i) * ((n - i) - 1) / 2 + j - i - 1;
            distance_matrix[index] = (*this->distance_metric_)(input_strings[i], input_strings[j]);
        }
    }

    return distance_matrix;
}
