//
// Created by Benjamin Tenmann on 19/11/2021.
//

#include <numeric>
#include <vector>
#include "metrics/Levenshtein.h"

double Levenshtein::forward(const std::string &a, const std::string &b) {
    size_t N = a.size();
    size_t M = b.size();

    if (!std::min(N, M)) return std::max(N, M);

    const char *c = &a.front();
    const char *d = &b.front();
    while (N > 0 && M > 0 && *c == *d) {
        c++;
        d++;
        N--;
        M--;
    }

    c = &a.back();
    d = &b.back();
    while (N > 0 && M > 0 && *c == *d) {
        c--;
        d--;
        N--;
        M--;
    }

    // below we set K to minimum and M to maximum
    size_t K = std::min(N, M);
    if (!K) return std::max(N, M);

    const std::string *e = &a;
    const std::string *f = &b;
    if (K == M) {
        M = N;
        e = &b;
        f = &a;
    }

    if (K == 1) {
        return (double) M - ((*f).find((*e).front()) != std::string::npos);
    }
    K++;
    M++;
    size_t half = K >> 1;
    std::vector<size_t> row(M);
    std::iota(row.begin(), row.end() - half, 0);

    // error code (?)
    if (row.empty()) return -1;

    return 0;
}
