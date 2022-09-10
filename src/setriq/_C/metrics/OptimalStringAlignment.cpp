//
// Created by Benjamin Tenmann on 10/06/2022.
//

#include "metrics/OptimalStringAlignment.h"

inline size_t min(const size_t& deletion, const size_t& insertion, const size_t& substitution) {
    if ((deletion < insertion) && (deletion < substitution))
        return deletion;
    else if ((insertion < deletion) && (insertion < substitution))
        return insertion;
    return substitution;
}

double metric::OptimalStringAlignment::forward(const std::string &a, const std::string &b) {
    const auto& len_a = a.size();
    const auto& len_b = b.size();

    // catch trivial cases
    if (!len_a) return (double) len_b;
    if (!len_b) return (double) len_a;

    auto&& score_matrix = uint_matrix_t(len_a + 1, uint_vector_t(len_b + 1, 0));
    for (auto i = 0ul; i <= len_a; i++) {
        for (auto j = 0ul; j <= len_b; j++) {
            if ((!i) || (!j)) {
                score_matrix[i][j] = (i + j);
                continue;
            }
            const auto& cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
            score_matrix[i][j] = min(
                    score_matrix[i - 1][j] + 1,
                    score_matrix[i][j - 1] + 1,
                    score_matrix[i - 1][j - 1] + cost
                    );
            if ((i > 1) && (j > 1) && (a[i - 1] == b[j]) && (a[i] == b[j - 1])) {
                const auto &transposition = score_matrix[i - 2][j - 2] + 1;
                if (transposition < score_matrix[i][j])
                    score_matrix[i][j] = transposition;
            }
        }
    }
    return (double) score_matrix[len_a][len_b];
}
