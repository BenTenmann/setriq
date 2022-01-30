//
// Created by Benjamin Tenmann on 19/11/2021.
//

#ifndef METRICS_TCRDIST_H
#define METRICS_TCRDIST_H

#include <string>
#include <vector>

#include "alignment/SubstitutionMatrix.h"

namespace metric {
    class TcrDist {
    private:
        SubstitutionMatrix substitution_matrix_;
        double gap_penalty_;
        char gap_symbol_;
        double distance_weight_;

    public:
        TcrDist() : substitution_matrix_{}, gap_penalty_{}, gap_symbol_{}, distance_weight_{} {};
        TcrDist(const double_matrix_t &scoring_matrix,
                const token_index_map_t &index,
                double gap_penalty,
                char gap_symbol,
                double weight);

        double forward(const std::string &, const std::string &) const;
    };
}

#endif //METRICS_TCRDIST_H
