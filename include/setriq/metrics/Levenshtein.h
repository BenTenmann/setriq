//
// Created by Benjamin Tenmann on 05/12/2021.
//

#ifndef SETRIQ_LEVENSHTEIN_H
#define SETRIQ_LEVENSHTEIN_H

#include <string>
#include "Metric.h"

namespace metric {
    class Levenshtein : public Metric {
        double extra_cost_ {0};

    public:
        Levenshtein() : extra_cost_ {} {};
        explicit Levenshtein(double xCost) : extra_cost_ {xCost} {};

        double forward(const std::string &, const std::string &) override;

    };
}

#endif //SETRIQ_LEVENSHTEIN_H
