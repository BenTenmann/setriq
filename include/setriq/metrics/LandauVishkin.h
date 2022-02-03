//
// Created by Benjamin Tenmann on 03/02/2022.
//

#ifndef SETRIQ_LANDAUVISHKIN_H
#define SETRIQ_LANDAUVISHKIN_H

#include "utils/type_defs.h"

namespace metric {
    class LandauVishkin {
    public:
        LandauVishkin() = default;

        double forward(const std::string&, const std::string&) const;
    };
}

#endif //SETRIQ_LANDAUVISHKIN_H
