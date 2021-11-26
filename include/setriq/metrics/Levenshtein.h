//
// Created by Benjamin Tenmann on 19/11/2021.
//

#ifndef METRICS_LEVENSHTEIN_H
#define METRICS_LEVENSHTEIN_H

#include <string>
#include "metrics/Metric.h"

class Levenshtein : public Metric {
public:
    double forward(const std::string &, const std::string&) override;
};


#endif //METRICS_LEVENSHTEIN_H
