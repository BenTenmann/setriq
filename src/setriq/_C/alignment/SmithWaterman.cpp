//
// Created by Benjamin Tenmann on 21/11/2021.
//


#include <utility>
#include "alignment/SmithWaterman.h"

SmithWaterman::SmithWaterman(SubstitutionMatrix matrix, double gapPen) {
    /**
     * Initialize a SmithWaterman object.
     *
     * @param matrix: a SubstitutionMatrix objet defining the substitution scores in the alignment
     * @param gapPen: the penalty for a gap in the alignment
     */
    this->substitutionMatrix = std::move(matrix);
    this->gapPenalty = gapPen;
}

double SmithWaterman::fillScoringMatrix(doubleMatrix &scoringMatrix,
                                        const std::string &a,
                                        const std::string &b) {
    /**
     * Fill the alignment scoring matrix and return the maximal alignment score between two input strings.
     *
     * @param scoringMatrix: the unfilled scoring matrix
     * @param a: an input string to be aligned
     * @param b: an input string to be aligned
     * @return the maximal alignment score between two sequences
     */
    size_t N = scoringMatrix.size();
    size_t M = scoringMatrix[0].size();

    double alignmentScore, kGapScore, lGapScore;
    double currentScore, maxScore {0};
    size_t kAxis = 0, lAxis = 1;

    for (size_t i = 1; i < N; i++) {
        for (size_t j = 1; j < M; j++) {
            alignmentScore = scoringMatrix[i - 1][j - 1] + this->substitutionMatrix(a[i - 1], b[j - 1]);
            kGapScore = this->calculateGapPenalty(scoringMatrix, i, j, kAxis);
            lGapScore = this->calculateGapPenalty(scoringMatrix, j, i, lAxis);

            currentScore = std::max(alignmentScore, std::max(kGapScore, lGapScore));
            if (currentScore > maxScore) maxScore = currentScore;
            scoringMatrix[i][j] = currentScore;
        }
    }
    return maxScore;
}

double SmithWaterman::calculateGapPenalty(const doubleMatrix &scoringMatrix,
                                          const size_t &maxGapLength,
                                          const size_t &index,
                                          const size_t &axis) const {
    /**
     * Calculate the gap score along a give axis. It has a lower bound of 0 and is linear (non-affine).
     *
     * @param scoringMatrix: the alignment scoring matrix
     * @param maxGapLength: the maximum possible gap length
     * @param index: the current index position (row or column)
     * @param axis: the axis over which to iterate
     * @return the maximum gap-score
     */
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

double SmithWaterman::computeBestAlignmentScore(const std::string &a, const std::string &b) {
    /**
     * Compute the maximal alignment score between two strings
     *
     * @param a: an input string to be aligned
     * @param b: an input string to be aligned
     * @return the maximal alignment score
     */
    size_t N = a.size();
    size_t M = b.size();

    doubleMatrix scoringMatrix (N + 1, doubleVector (M + 1, 0));
    double score {this->fillScoringMatrix(scoringMatrix, a, b)};

    return score;
}

double SmithWaterman::forward(const std::string &a, const std::string &b) {
    /**
     * Wrapper for computing the maximal alignment score
     *
     * @param a: an input string to be aligned
     * @param b: an input string to be aligned
     * @return the maximal alignment score
     */
    double bestScore {this->computeBestAlignmentScore(a, b)};
    return bestScore;
}

double SmithWaterman::identityScore(const std::string &inputString) {
    /**
     * Compute the maximal alignment score for an input string with itself. It is a trivial case of SW, as the maximal
     * score will always be at the final diagonal position of the scoring matrix. Thus it reduces down to a cumulative
     * sum of the substitution scores of the string characters with themselves.
     *
     * @param inputString: an input string to be aligned with itself
     * @return the maximal self-alignment score for an input string
     */
    size_t N {inputString.size()};

    double score {0};
    for (size_t i = 0; i < N; i++) {
        score += this->substitutionMatrix(inputString[i], inputString[i]);
    }
    return score;
}
