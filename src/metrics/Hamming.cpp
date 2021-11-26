//
// Created by Benjamin Tenmann on 19/11/2021.
//

#include <algorithm>
#include <cmath>
#include "metrics/Hamming.h"

double Hamming::forward(const std::string &a, const std::string &b) {
    unsigned int dist = 0;
    size_t N = a.size();
    size_t M = b.size();

    double leftOver = std::abs(int(N - M));
    size_t minStringLength {std::min(N, M)};
    for (size_t index = 0; index < minStringLength; index++) {
        dist += a[index] != b[index];
    }
    return dist + leftOver;
}
