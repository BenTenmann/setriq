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
     * @return a flat (N * (N - 1)) / 2 vector of doubles
     */
    const auto& N = input_strings.size();
    double_vector_t distance_matrix(N * (N - 1) / 2);

    #pragma omp parallel for default(none) shared(input_strings, distance_matrix, N)
    for (size_t i = 0; i < N - 1; i++) {
        for (size_t j = i + 1; j < N; j++) {
            size_t index = (N * (N - 1)) / 2 - (N - i) * ((N - i) - 1) / 2 + j - i - 1;
            distance_matrix[index] = (*this->distance_metric_)(input_strings[i], input_strings[j]);
        }
    }

    return distance_matrix;
}
