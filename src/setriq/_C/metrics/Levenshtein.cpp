//
// Created by Benjamin Tenmann on 05/12/2021.
//

#include <cstring>
#include <numeric>
#include <string>
#include <vector>
#include "metrics/Levenshtein.h"

double metric::Levenshtein::forward(const std::string &a, const std::string &b) {
    size_t lengthOfA {a.size()};
    size_t lengthOfB {b.size()};
    size_t halfOfLengthA;

    const char* ptrToA = &a.front();
    const char* ptrToB = &b.front();

    // catch the trivial cases
    if (a.empty()) return lengthOfB;
    if (b.empty()) return lengthOfA;

    // grind down common prefix
    while (lengthOfA > 0 && lengthOfB > 0 && (*ptrToA) == (*ptrToB)) {
        lengthOfA--;
        lengthOfB--;
        ptrToA++;
        ptrToB++;
    }

    // grind down common suffix
    while (lengthOfA > 0 && lengthOfB > 0 && ptrToA[lengthOfA - 1] == ptrToB[lengthOfB - 1]) {
        lengthOfA--;
        lengthOfB--;
    }

    // again, catch trivial cases
    if (lengthOfA == 0) return lengthOfB;
    if (lengthOfB == 0) return lengthOfA;

    if (lengthOfA > lengthOfB) {  // enforce that b is the longer string
        size_t temporaryLengthStore = lengthOfA;
        const char *temporaryPtrStore = ptrToA;
        lengthOfA = lengthOfB;
        lengthOfB = temporaryLengthStore;
        ptrToA = ptrToB;
        ptrToB = temporaryPtrStore;
    }

    if (lengthOfA == 1) {
        if (this->extraCost > 0)
            return (double) lengthOfB + 1 - this->extraCost * (memchr(ptrToB, *ptrToA, lengthOfB) != nullptr);
        else
            return (double) lengthOfB - (memchr(ptrToB, *ptrToA, lengthOfB) != nullptr);
    }
    lengthOfA++;
    lengthOfB++;
    halfOfLengthA = lengthOfA >> 1;

    // first row initialization
    std::vector<size_t> row (lengthOfB);
    std::iota(row.begin(), row.end() - (this->extraCost > 0 ? 0 : halfOfLengthA), 0);

    size_t rowIndex;
    size_t *end = &row.back() - lengthOfB - 1;

    if (this->extraCost > 0) {
        for (rowIndex = 1; rowIndex < lengthOfA; rowIndex++) {
            size_t *ptrToRowElement = &row[1];
            const char currentCharFromA = ptrToA[rowIndex - 1];
            const char *ptrToCurrentCharFromB = ptrToB;

            size_t rowIndexCopy1 = rowIndex;
            size_t rowIndexCopy2 = rowIndex;
            while (ptrToRowElement <= end) {
                if (currentCharFromA == *(ptrToCurrentCharFromB++))
                    rowIndexCopy2 = --rowIndexCopy1;
                else
                    rowIndexCopy2++;

                rowIndexCopy1 = *ptrToRowElement;
                rowIndexCopy1++;

                if (rowIndexCopy2 > rowIndexCopy1)
                    rowIndexCopy2 = rowIndexCopy1;

                *(ptrToRowElement++) = rowIndexCopy2;
            }
        }
    }

    else {
        /*
         *
         * in this case we don't have to scan two corner triangles (of size len1/2)
         * in the matrix because no best path can go throught them. note this
         * breaks when len1 == len2 == 2 so the memchr() special case above is
         * necessary
         *
         */
        row[0] = lengthOfA - halfOfLengthA - 1;
        for (rowIndex = 1; rowIndex < lengthOfA; rowIndex++) {
            size_t *ptrToRowElement;
            const char currentCharFromA = ptrToA[rowIndex - 1];
            const char *ptrToCurrentCharFromB;

            size_t rowIndexCopy1, rowIndexCopy2;
            /* skip the upper triangle */
            if (rowIndex >= lengthOfA - halfOfLengthA) {
                size_t offset = rowIndex - (lengthOfA - halfOfLengthA);
                size_t rowIndexCopy3;

                ptrToCurrentCharFromB = ptrToB + offset;
                ptrToRowElement = &row[offset];

                rowIndexCopy3 = *(ptrToRowElement++) + (currentCharFromA != *(ptrToCurrentCharFromB++));
                rowIndexCopy2 = *ptrToRowElement;
                rowIndexCopy2++;
                rowIndexCopy1 = rowIndexCopy2;
                if (rowIndexCopy2 > rowIndexCopy3)
                    rowIndexCopy2 = rowIndexCopy3;
                *(ptrToRowElement++) = rowIndexCopy2;
            }
            else {
                ptrToRowElement = &row[1];
                ptrToCurrentCharFromB = ptrToB;
                rowIndexCopy1 = rowIndexCopy2 = rowIndex;
            }

            /* skip the lower triangle */
            if (rowIndex <= halfOfLengthA + 1)
                end = &row[lengthOfB + rowIndex - halfOfLengthA - 2];

            /* main */
            while (ptrToRowElement <= end) {
                size_t rowIndexCopy3 = --rowIndexCopy2 + (currentCharFromA != *(ptrToCurrentCharFromB++));
                rowIndexCopy1++;

                if (rowIndexCopy1 > rowIndexCopy3)
                    rowIndexCopy1 = rowIndexCopy3;

                rowIndexCopy2 = *ptrToRowElement;
                rowIndexCopy2++;

                if (rowIndexCopy1 > rowIndexCopy2)
                    rowIndexCopy1 = rowIndexCopy2;
                *(ptrToRowElement++) = rowIndexCopy1;
            }

            /* lower triangle sentinel */
            if (rowIndex <= halfOfLengthA) {
                size_t rowIndexCopy3 = --rowIndexCopy1 + (currentCharFromA != (*ptrToCurrentCharFromB));
                rowIndexCopy2++;
                if (rowIndexCopy2 > rowIndexCopy3)
                    rowIndexCopy2 = rowIndexCopy3;
                *ptrToRowElement = rowIndexCopy2;
            }
        }
    }

    return (double) *end;
}
