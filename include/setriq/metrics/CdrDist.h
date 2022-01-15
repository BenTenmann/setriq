//
// Created by Benjamin Tenmann on 20/11/2021.
//

#ifndef METRICS_CDRDIST_H
#define METRICS_CDRDIST_H

#include "Metric.h"
#include "alignment/SmithWaterman.h"

namespace metric {

    class CdrDist : public Metric {
    private:
        SmithWaterman algorithm_;

    public:
        CdrDist() : algorithm_{} {};
        explicit CdrDist(const double_matrix_t&, const token_index_map_t&, const double&, const double&);

        double forward(const std::string &, const std::string &) override;
    };

}

#endif //METRICS_CDRDIST_H
