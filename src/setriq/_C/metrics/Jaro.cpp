//
// Created by Benjamin Tenmann on 05/03/2022.
//

#include <cmath>
#include "metrics/Jaro.h"

#define either_zero(x, y) (x == 0) || (y == 0)
#define max(x, y) x > y ? x : y
#define min(x, y) x > y ? y : x


void collapse_into_match_str(const std::string& sequence, const std::vector<size_t>& matches_idx, char* match_str) {
    auto&& j = 0ul;
    for (const auto& idx : matches_idx) {
        if (idx){
            match_str[j] = sequence[idx - 1];
            j++;
        }
    }
}

double metric::Jaro::forward(const std::string &a, const std::string &b) const {
    /*!
     * Compute the Jaro distance between two input strings.
     * Adapted from https://github.com/markvanderloo/stringdist/blob/master/pkg/src/jaro.c
     *
     * @param a: an input string to be compared
     * @param b: an input string to be compared
     */
    const auto& s_i = a.size();
    const auto& s_j = b.size();
    if (either_zero(s_i, s_j))
        // if one of the strings is of length 0 and the other isn't, then the distance is maximal (1)
        // if both are length 0, then the distance is minimal, i.e. 0
        return (double) ((s_i > 0) || (s_j > 0));

    const auto& max_len = s_i > s_j ? s_i : s_j;
    const auto& max_match_distance = (int) std::floor(max_len / 2) - 1;
    if (max_match_distance < 0)
        // catch the case when both strings are of length == 1
        return a[0] == b[0] ? 0.0 : 1.0;

    auto&& matches_s_i = std::vector<size_t>(s_i, 0);
    auto&& matches_s_j = std::vector<size_t>(s_j, 0);

    auto&& n_matches = 0ul;
    for (auto i = 0; i < s_i; i++) {
        const auto& left = max((i - max_match_distance), 0);
        const auto& right = min((i + max_match_distance) + 1, s_j);
        // can we collapse this in some way?
        for (auto j = left; j < right; j++) {
            if ((a[i] == b[j]) && (matches_s_j[j] == 0)) {
                n_matches++;
                matches_s_i[i] = i + 1;
                matches_s_j[j] = j + 1;
                break;
            }
        }
    }
    if (n_matches == 0)
        return 1.0;

    char *match_str_i = new char[n_matches];
    char *match_str_j = new char[n_matches];

    collapse_into_match_str(a, matches_s_i, match_str_i);
    collapse_into_match_str(b, matches_s_j, match_str_j);

    auto&& t = 0.0;
    for (auto k = 0ul; k < n_matches; k++) {
        if (match_str_i[k] != match_str_j[k])
            t += 0.5;
    }
    delete []match_str_i;
    delete []match_str_j;

    const auto& m = (double) n_matches;
    // allow arbitrary weighting
    return 1 - (this->weights_[0] * (m / s_i) + this->weights_[1] * (m / s_j) + this->weights_[2] * ((m - t) / m));
}
