#include <__bit_reference>
//
// Created by Benjamin Tenmann on 19/11/2021.
//

#ifndef METRICS_TCRDIST_H
#define METRICS_TCRDIST_H

#include <string>
#include <vector>
#include "metrics/Metric.h"
#include "alignment/SubstitutionMatrix.h"

class TcrDist : public Metric {
private:
    __unused SubstitutionMatrix substitutionMatrix;

public:
    double forward(const std::string&, const std::string&) override;
};


#endif //METRICS_TCRDIST_H
