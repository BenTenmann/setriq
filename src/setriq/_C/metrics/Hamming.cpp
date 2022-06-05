//
// Created by Benjamin Tenmann on 20/02/2022.
//

#include "metrics/Hamming.h"

double metric::Hamming::forward(const std::string &a, const std::string &b) const {
    /*!
     * Compute the Hamming distance between two input strings.
     *
     * @param a: an input string to be compared
     * @param b: an input string to be compared
     */
    auto&& distance = 0.;
    for (auto i = 0ul; i < a.size(); i++) {
        if (a[i] != b[i])
            distance += this->mismatch_score_;
    }
    return distance;
}
