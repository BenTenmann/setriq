//
// Created by Benjamin Tenmann on 19/11/2021.
//

#ifndef METRICS_TCRDIST_H
#define METRICS_TCRDIST_H

#include <string>
#include <vector>
#include "metrics/Metric.h"
#include "alignment/SubstitutionMatrix.h"

namespace metric {
    class TcrDist : public Metric {
    private:
        SubstitutionMatrix substitutionMatrix;
        double gapPenalty;
        char gapSymbol;
        double distanceWeight;

    public:
        TcrDist() : substitutionMatrix(), gapPenalty(), gapSymbol(), distanceWeight() {};
        TcrDist(const doubleMatrix &subMat, const stringIndexMap &, double gapPen, char gapSym, double weight);

        double forward(const std::string &, const std::string &) override;
    };
}

#endif //METRICS_TCRDIST_H
