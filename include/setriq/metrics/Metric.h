//
// Created by Benjamin Tenmann on 19/11/2021.
//

#ifndef METRICS_METRIC_H
#define METRICS_METRIC_H

#include <string>

class Metric {
public:
    Metric() = default;

    virtual double forward(const std::string&, const std::string&) = 0;
    double operator () (const std::string &a, const std::string&b) {
        double out {this->forward(a, b)};

        return out;
    };
};


#endif //METRICS_METRIC_H
