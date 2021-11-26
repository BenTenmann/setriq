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
        SmithWaterman algorithm;

    public:
        CdrDist() : algorithm() {};
        explicit CdrDist(const doubleMatrix&, const stringIndexMap&, double = 1.0, size_t = 1000);

        double forward(const std::string &, const std::string &) override;
    };

}

#endif //METRICS_CDRDIST_H
