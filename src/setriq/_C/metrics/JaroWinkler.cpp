//
// Created by Benjamin Tenmann on 21/02/2022.
//

#include "metrics/JaroWinkler.h"

size_t min_sequence_len(const std::string& a, const std::string& b) {
    const auto& length_a = a.size();
    const auto& length_b = b.size();
    return length_a < length_b ? length_a : length_b;
}

double metric::JaroWinkler::forward(const std::string &a, const std::string &b) const {
    /*!
     * Compute the Jaro-Winkler distance between two input strings.
     *
     * @param a: an input string to be compared
     * @param b: an input string to be compared
     */
    const auto& jaro_distance = this->jaro_.forward(a, b);
    const auto& min_length = min_sequence_len(a, b);

    auto&& l = 0ul;
    while ((a[l] == b[l]) && (l < min_length) && (l < this->max_l_))
        l++;
    return jaro_distance * (1 - l * this->p_);
}
