//
// Created by Benjamin Tenmann on 21/11/2021.
//

#ifndef METRICS_SMITHWATERMAN_H
#define METRICS_SMITHWATERMAN_H

#include "SubstitutionMatrix.h"
#include "utils/typeDefs.h"

class SmithWaterman {
private:
    SubstitutionMatrix substitution_matrix_;
    double gap_opening_penalty_;
    double gap_extension_penalty_;

    // scoring matrix creation
    double calculate_gap_penalty_(const double_matrix_t& scoring_matrix,
                                  const size_t& max_gap_length,
                                  const size_t& index,
                                  const size_t& axis) const;
    double fill_scoring_matrix_(double_matrix_t& scoring_matrix, const std::string&a, const std::string&b);
    double compute_best_alignment_score_(const std::string&a, const std::string&b);

public:
    SmithWaterman() : substitution_matrix_{}, gap_opening_penalty_{}, gap_extension_penalty_{} {};
    SmithWaterman(SubstitutionMatrix, const double&, const double&);

    double forward(const std::string&, const std::string&);
    double operator() (const std::string& a, const std::string& b) {
        double out {this->forward(a, b)};

        return out;
    };

    // special case: identity score
    double identity_score(const std::string&);
};


#endif //METRICS_SMITHWATERMAN_H
