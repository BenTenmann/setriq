//
// Created by Benjamin Tenmann on 20/11/2021.
//

#ifndef METRICS_CDRDIST_H
#define METRICS_CDRDIST_H

#include "alignment/SmithWaterman.h"

namespace metric {

    class CdrDist {
    private:
        SmithWaterman algorithm_;

    public:
        CdrDist() : algorithm_{} {};
        explicit CdrDist(const double_matrix_t&, const token_index_map_t&, const double&, const double&);

        double forward(const std::string &, const std::string &) const;
    };

}

#endif //METRICS_CDRDIST_H
