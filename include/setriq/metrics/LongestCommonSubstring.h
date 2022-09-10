//
// Created by Benjamin Tenmann on 10/06/2022.
//

#ifndef SETRIQ_LONGESTCOMMONSUBSTRING_H
#define SETRIQ_LONGESTCOMMONSUBSTRING_H

#include "utils/type_defs.h"

namespace metric {
    class LongestCommonSubstring {
    public:
        LongestCommonSubstring() = default;

        double forward(const std::string&, const std::string&);
    };
}

#endif //SETRIQ_LONGESTCOMMONSUBSTRING_H
