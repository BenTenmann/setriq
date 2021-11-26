//
// Created by Benjamin Tenmann on 19/11/2021.
//

#ifndef METRICS_METRIC_H
#define METRICS_METRIC_H

#include <string>
#include "utils/LRUCache.h"

class Metric {
public:
    LRUCache <std::string, double> cache {1000};
    Metric() = default;

    virtual double forward(const std::string&, const std::string&) = 0;
    double operator () (const std::string &a, const std::string&b) {
        if (this->cache.exists(a + b)) return this->cache.get(a + b);
        if (this->cache.exists(b + a)) return this->cache.get(b + a);

        double out {this->forward(a, b)};
        this->cache.put(a + b, out);
        return out;
    };
};


#endif //METRICS_METRIC_H
