//
// Created by Benjamin Tenmann on 21/02/2022.
//

#ifndef SETRIQ_JAROWINKLER_H
#define SETRIQ_JAROWINKLER_H

#include "utils/type_defs.h"
#include "metrics/Jaro.h"

namespace metric {
    class JaroWinkler {
    private:
        double p_ = 0.;
        size_t max_l_ = 4;
        Jaro jaro_{};

    public:
        JaroWinkler() = default;

        explicit JaroWinkler(const double &p, const size_t &max_l, Jaro jaro) : p_{p}, max_l_{max_l}, jaro_{jaro} {};

        double forward(const std::string &a, const std::string &b) const;
    };
}

#endif //SETRIQ_JAROWINKLER_H
