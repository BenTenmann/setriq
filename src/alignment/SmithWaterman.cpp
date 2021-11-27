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

doubleMatrix SmithWaterman::initialiseScoringMatrix(const size_t &N, const size_t &M) {
    doubleMatrix scoringMatrix;
    scoringMatrix.reserve(N + 1);
    for (size_t i = 0; i < (N + 1); i++) {
        doubleVector row(M + 1, 0);
        scoringMatrix.push_back(row);
    }
    return scoringMatrix;
}

void SmithWaterman::fillScoringMatrix(doubleMatrix &scoringMatrix,
                                      const std::string &a,
                                      const std::string &b) {
    size_t N = scoringMatrix.size();
    size_t M = scoringMatrix[0].size();

    double alignmentScore, kGapScore, lGapScore;
    size_t kAxis = 0, lAxis = 1;

    for (size_t i = 1; i < N; i++) {
        for (size_t j = 1; j < M; j++) {
            alignmentScore = scoringMatrix[i - 1][j - 1] + this->substitutionMatrix(a.substr(i - 1, 1),
                                                                                    b.substr(j - 1, 1));
            kGapScore = this->calculateGapPenalty(scoringMatrix, i, j, kAxis);
            lGapScore = this->calculateGapPenalty(scoringMatrix, j, i, lAxis);

            scoringMatrix[i][j] = std::max(alignmentScore, std::max(kGapScore, lGapScore));
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

doubleMatrix SmithWaterman::createScoringMatrix(const std::string &a, const std::string &b) {
    size_t N = a.size();
    size_t M = b.size();

    doubleMatrix scoringMatrix = SmithWaterman::initialiseScoringMatrix(N, M);
    this->fillScoringMatrix(scoringMatrix, a, b);

    return scoringMatrix;
}

double SmithWaterman::forward(const std::string &a, const std::string &b) {
    doubleMatrix scoringMatrix = this->createScoringMatrix(a, b);

    IndexTuple argMax = SmithWaterman::argMax(scoringMatrix);

    return SmithWaterman::tracebackScore(scoringMatrix, argMax);
}

IndexTuple SmithWaterman::argMax(const doubleMatrix &scoringMatrix) {
    size_t N = scoringMatrix.size();
    size_t M = scoringMatrix[0].size();

    size_t iMax = 0, jMax = 0;
    double currentScore, maxScore;

    maxScore = 0;
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            currentScore = scoringMatrix[i][j];
            if (currentScore > maxScore) {
                maxScore = currentScore;
                iMax = i;
                jMax = j;
            }
        }
    }
    return {iMax, jMax};
}

double SmithWaterman::tracebackScore(const doubleMatrix &scoringMatrix, const IndexTuple &argMax) {
    double totalScore {scoringMatrix[argMax.i][argMax.j]};

    double currentScore {totalScore};
    IndexTuple currentPosition {argMax.i, argMax.j};
    while (currentScore > 0) {
        currentPosition = SmithWaterman::tracebackStep(scoringMatrix, currentPosition);
        currentScore = scoringMatrix[currentPosition.i][currentPosition.j];

        totalScore += currentScore;
    }

    return totalScore;
}

IndexTuple SmithWaterman::tracebackStep(const doubleMatrix &scoringMatrix, const IndexTuple &currentPosition) {
    size_t iMax {currentPosition.i}, jMax {currentPosition.j};
    double across {scoringMatrix[iMax][jMax - 1]};
    double diagonal {scoringMatrix[iMax - 1][jMax - 1]};
    double up {scoringMatrix[iMax - 1][jMax]};

    if (across > diagonal && across > up) jMax--;

    else if (up > diagonal && up > across) iMax--;

    else {
        iMax--;
        jMax--;
    }

    return {iMax, jMax};
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
