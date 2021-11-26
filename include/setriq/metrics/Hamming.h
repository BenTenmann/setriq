//
// Created by Benjamin Tenmann on 19/11/2021.
//

#ifndef METRICS_HAMMING_H
#define METRICS_HAMMING_H

#include <string>
#include "metrics/Metric.h"

class Hamming : public Metric {
public:
    double forward(const std::string&, const std::string&) override;
};


#endif //METRICS_HAMMING_H
