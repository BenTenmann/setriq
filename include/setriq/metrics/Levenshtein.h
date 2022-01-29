//
// Created by Benjamin Tenmann on 05/12/2021.
//

#ifndef SETRIQ_LEVENSHTEIN_H
#define SETRIQ_LEVENSHTEIN_H

#include <string>

namespace metric {
    class Levenshtein {
        double extra_cost_ = 0;

    public:
        Levenshtein() : extra_cost_{} {};
        explicit Levenshtein(double x_cost) : extra_cost_ {x_cost} {};

        double forward(const std::string &, const std::string &) const;

    };
}

#endif //SETRIQ_LEVENSHTEIN_H
