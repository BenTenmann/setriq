//
// Created by Benjamin Tenmann on 21/11/2021.
//

#ifndef METRICS_SMITHWATERMAN_H
#define METRICS_SMITHWATERMAN_H

#include "SubstitutionMatrix.h"
#include "utils/typeDefs.h"

class SmithWaterman {
private:
    SubstitutionMatrix substitutionMatrix;
    double gapPenalty;

    // scoring matrix creation
    double calculateGapPenalty(const doubleMatrix&, const size_t&, const size_t&, const size_t&) const;
    double fillScoringMatrix(doubleMatrix&, const std::string&, const std::string&);
    double computeBestAlignmentScore(const std::string&, const std::string&);

public:
    SmithWaterman() : substitutionMatrix(), gapPenalty{0} {};
    SmithWaterman(SubstitutionMatrix, double);

    double forward(const std::string&, const std::string&);
    double operator() (const std::string& a, const std::string& b) {
        double out {this->forward(a, b)};

        return out;
    };

    // special case: identity score
    double identityScore(const std::string&);
};


#endif //METRICS_SMITHWATERMAN_H
