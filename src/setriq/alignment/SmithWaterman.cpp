//
// Created by Benjamin Tenmann on 21/11/2021.
//

#include <utility>
#include "alignment/SmithWaterman.h"

SmithWaterman::SmithWaterman(SubstitutionMatrix matrix, double gapPen, size_t cacheSize) {
    this->cache = LRUCache<std::string, double> (cacheSize);
    this->substitutionMatrix = std::move(matrix);
    this->gapPenalty = gapPen;
}

void SmithWaterman::fillScoringMatrix(doubleMatrix &scoringMatrix,
                                      const std::string &a,
                                      const std::string &b,
                                      IndexTuple &argMax) {
    size_t N = scoringMatrix.size();
    size_t M = scoringMatrix[0].size();

    double alignmentScore, kGapScore, lGapScore;
    double currentScore, maxScore {0};
    size_t kAxis = 0, lAxis = 1;

    for (size_t i = 1; i < N; i++) {
        for (size_t j = 1; j < M; j++) {
            alignmentScore = scoringMatrix[i - 1][j - 1] + this->substitutionMatrix(a.substr(i - 1, 1),
                                                                                    b.substr(j - 1, 1));
            kGapScore = this->calculateGapPenalty(scoringMatrix, i, j, kAxis);
            lGapScore = this->calculateGapPenalty(scoringMatrix, j, i, lAxis);

            currentScore = std::max(alignmentScore, std::max(kGapScore, lGapScore));
            if (currentScore > maxScore) {
                maxScore = currentScore;
                argMax.i = i;
                argMax.j = j;
            }
            scoringMatrix[i][j] = currentScore;
        }
    }
}

double SmithWaterman::calculateGapPenalty(const doubleMatrix &scoringMatrix,
                                          const size_t &maxGapLength,
                                          const size_t &index,
                                          const size_t &axis) const {
    double maxScore = 0;

    double elem, score;
    size_t k;
    for (size_t i = 1; i < (maxGapLength + 1); i++) {
        k = maxGapLength - i;
        elem = axis ? scoringMatrix[index][k] : scoringMatrix[k][index];
        score = elem - (i * this->gapPenalty);

        if (score > maxScore) maxScore = score;
    }
    return maxScore;
}

doubleMatrix SmithWaterman::createScoringMatrix(const std::string &a, const std::string &b, IndexTuple &argMax) {
    size_t N = a.size();
    size_t M = b.size();

    doubleMatrix scoringMatrix (N + 1, doubleVector (M + 1, 0));
    this->fillScoringMatrix(scoringMatrix, a, b, argMax);

    return scoringMatrix;
}

double SmithWaterman::forward(const std::string &a, const std::string &b) {
    IndexTuple argMax {0, 0};
    doubleMatrix scoringMatrix {this->createScoringMatrix(a, b, argMax)};

    return SmithWaterman::tracebackScore(scoringMatrix, argMax);
}

double SmithWaterman::tracebackScore(const doubleMatrix &scoringMatrix, IndexTuple &currentPosition) {
    double totalScore {scoringMatrix[currentPosition.i][currentPosition.j]};

    double currentScore {totalScore};
    while (currentScore > 0) {
        SmithWaterman::tracebackStep(scoringMatrix, currentPosition.i, currentPosition.j);
        currentScore = scoringMatrix[currentPosition.i][currentPosition.j];

        totalScore += currentScore;
    }

    return totalScore;
}

void SmithWaterman::tracebackStep(const doubleMatrix &scoringMatrix, size_t &i, size_t &j) {
    double across {scoringMatrix[i][j - 1]};
    double diagonal {scoringMatrix[i - 1][j - 1]};
    double up {scoringMatrix[i - 1][j]};

    if (across > diagonal && across > up) j--;

    else if (up > diagonal && up > across) i--;

    else {
        i--;
        j--;
    }
}

double SmithWaterman::identityScore(const std::string &inputString) {
    size_t N {inputString.size()};

    double score {0};
    for (size_t i = 0; i < N; i++) {
        double elem = this->substitutionMatrix(inputString.substr(i, 1), inputString.substr(i, 1));
        score += (double) (N - i) * elem;
    }
    return score;
}
