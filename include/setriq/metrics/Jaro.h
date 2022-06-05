//
// Created by Benjamin Tenmann on 05/03/2022.
//

#ifndef SETRIQ_JARO_H
#define SETRIQ_JARO_H

#include <array>
#include "utils/type_defs.h"

typedef std::array<double, 3> jaro_weighting_t;

namespace metric {
    class Jaro {
    private:
        jaro_weighting_t weights_ = {1. / 3, 1. / 3, 1. / 3};

    public:
        Jaro() = default;

        explicit Jaro(jaro_weighting_t weights) : weights_(weights) {};

        double forward(const std::string &a, const std::string &b) const;
    };
}

#endif //SETRIQ_JARO_H
