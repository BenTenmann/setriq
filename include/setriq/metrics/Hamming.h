//
// Created by Benjamin Tenmann on 20/02/2022.
//

#ifndef SETRIQ_HAMMING_H
#define SETRIQ_HAMMING_H

#include "utils/type_defs.h"

namespace metric {
    class Hamming {
    private:
        double mismatch_score_{};

    public:
        explicit Hamming(const double &mismatch_score) : mismatch_score_{mismatch_score} {};

        double forward(const std::string &a, const std::string &b) const;
    };
}

#endif //SETRIQ_HAMMING_H
