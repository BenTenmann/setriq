//
// Created by Benjamin Tenmann on 05/12/2021.
//

#ifndef SETRIQ_LEVENSHTEIN_H
#define SETRIQ_LEVENSHTEIN_H

#include <string>
#include "Metric.h"

namespace metric {
    class Levenshtein : public Metric {
        double extraCost = 0;

    public:
        Levenshtein() : extraCost() {};
        explicit Levenshtein(double xCost) : extraCost {xCost} {};

        double forward(const std::string &, const std::string &);

    };
}

#endif //SETRIQ_LEVENSHTEIN_H
