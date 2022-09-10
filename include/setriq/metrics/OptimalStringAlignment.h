//
// Created by Benjamin Tenmann on 10/06/2022.
//

#ifndef SETRIQ_OPTIMALSTRINGALIGNMENT_H
#define SETRIQ_OPTIMALSTRINGALIGNMENT_H

#include "utils/type_defs.h"

namespace metric {
    class OptimalStringAlignment {
    public:
        OptimalStringAlignment() = default;

        double forward(const std::string&, const std::string&);
    };
}

#endif //SETRIQ_OPTIMALSTRINGALIGNMENT_H
