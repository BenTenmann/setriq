//
// Created by Benjamin Tenmann on 05/12/2021.
//

#include <cstring>
#include <numeric>
#include <vector>

#include "metrics/Levenshtein.h"

double metric::Levenshtein::forward(const std::string &a, const std::string &b) const {
    /**
     * Compute the Levenshtein distance between two input strings. This is a C++ refactor of the `python-Levenshtein`
     * implementation (https://github.com/ztane/python-Levenshtein).
     *
     * @param a: an input string to be compared
     * @param b: an input string to be compared
     * @return the Levensthein distance between the two input strings
     */
    size_t length_of_a {a.size()};
    size_t length_of_b {b.size()};
    size_t half_of_length_a;

    const char* ptr_to_a = &a.front();
    const char* ptr_to_b = &b.front();

    // catch the trivial cases
    if (a.empty()) return length_of_b;
    if (b.empty()) return length_of_a;

    // grind down common prefix
    while (length_of_a > 0 && length_of_b > 0 && (*ptr_to_a) == (*ptr_to_b)) {
        length_of_a--;
        length_of_b--;
        ptr_to_a++;
        ptr_to_b++;
    }

    // grind down common suffix
    while (length_of_a > 0 && length_of_b > 0 && ptr_to_a[length_of_a - 1] == ptr_to_b[length_of_b - 1]) {
        length_of_a--;
        length_of_b--;
    }

    // again, catch trivial cases
    if (length_of_a == 0) return length_of_b;
    if (length_of_b == 0) return length_of_a;

    if (length_of_a > length_of_b) {  // enforce that b is the longer string
        size_t temporary_length_store = length_of_a;
        const char *temporary_ptr_store = ptr_to_a;
        length_of_a = length_of_b;
        length_of_b = temporary_length_store;
        ptr_to_a = ptr_to_b;
        ptr_to_b = temporary_ptr_store;
    }

    if (length_of_a == 1) {
        if (this->extra_cost_ > 0)
            return (double) length_of_b + 1 - this->extra_cost_ * (memchr(ptr_to_b, *ptr_to_a, length_of_b) != nullptr);
        else
            return (double) length_of_b - (memchr(ptr_to_b, *ptr_to_a, length_of_b) != nullptr);
    }
    length_of_a++;
    length_of_b++;
    half_of_length_a = length_of_a >> 1;

    // first row initialization
    std::vector<size_t> row (length_of_b);
    std::iota(row.begin(), row.end() - (this->extra_cost_ > 0 ? 0 : half_of_length_a), 0);

    size_t row_index;
    size_t *end = &row.back() - length_of_b - 1;

    if (this->extra_cost_ > 0) {
        for (row_index = 1; row_index < length_of_a; row_index++) {
            size_t *ptr_to_row_element = &row[1];
            const char current_char_from_a = ptr_to_a[row_index - 1];
            const char *ptr_to_current_char_from_b = ptr_to_b;

            size_t row_index_copy_1 = row_index;
            size_t row_index_copy_2 = row_index;
            while (ptr_to_row_element <= end) {
                if (current_char_from_a == *(ptr_to_current_char_from_b++))
                    row_index_copy_2 = --row_index_copy_1;
                else
                    row_index_copy_2++;

                row_index_copy_1 = *ptr_to_row_element;
                row_index_copy_1++;

                if (row_index_copy_2 > row_index_copy_1)
                    row_index_copy_2 = row_index_copy_1;

                *(ptr_to_row_element++) = row_index_copy_2;
            }
        }
    }

    else {
        /*
         *
         * in this case we don't have to scan two corner triangles (of size len1/2)
         * in the matrix because no best path can go through them. note this
         * breaks when length_of_a == length_of_b == 2 so the `memchr()` special case above is
         * necessary
         *
         */
        row[0] = length_of_a - half_of_length_a - 1;
        for (row_index = 1; row_index < length_of_a; row_index++) {
            size_t *ptr_to_row_element;
            const char current_char_from_a = ptr_to_a[row_index - 1];
            const char *ptr_to_current_char_from_b;

            size_t row_index_copy_1, row_index_copy_2;
            /* skip the upper triangle */
            if (row_index >= length_of_a - half_of_length_a) {
                size_t offset = row_index - (length_of_a - half_of_length_a);
                size_t row_index_copy_3;

                ptr_to_current_char_from_b = ptr_to_b + offset;
                ptr_to_row_element = &row[offset];

                row_index_copy_3 = *(ptr_to_row_element++) + (current_char_from_a != *(ptr_to_current_char_from_b++));
                row_index_copy_2 = *ptr_to_row_element;
                row_index_copy_2++;

                row_index_copy_1 = row_index_copy_2;
                if (row_index_copy_2 > row_index_copy_3)
                    row_index_copy_2 = row_index_copy_3;

                *(ptr_to_row_element++) = row_index_copy_2;
            }
            else {
                ptr_to_row_element = &row[1];
                ptr_to_current_char_from_b = ptr_to_b;
                row_index_copy_1 = row_index_copy_2 = row_index;
            }

            /* skip the lower triangle */
            if (row_index <= half_of_length_a + 1)
                end = &row[length_of_b + row_index - half_of_length_a - 2];

            /* main */
            while (ptr_to_row_element <= end) {
                size_t row_index_copy_3 = --row_index_copy_2 + (current_char_from_a != *(ptr_to_current_char_from_b++));
                row_index_copy_1++;

                if (row_index_copy_1 > row_index_copy_3)
                    row_index_copy_1 = row_index_copy_3;

                row_index_copy_2 = *ptr_to_row_element;
                row_index_copy_2++;

                if (row_index_copy_1 > row_index_copy_2)
                    row_index_copy_1 = row_index_copy_2;
                *(ptr_to_row_element++) = row_index_copy_1;
            }

            /* lower triangle sentinel */
            if (row_index <= half_of_length_a) {
                size_t row_index_copy_3 = --row_index_copy_1 + (current_char_from_a != (*ptr_to_current_char_from_b));

                row_index_copy_2++;
                if (row_index_copy_2 > row_index_copy_3)
                    row_index_copy_2 = row_index_copy_3;

                *ptr_to_row_element = row_index_copy_2;
            }
        }
    }

    return (double) *end;
}
