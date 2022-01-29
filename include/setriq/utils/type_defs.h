//
// Created by Benjamin Tenmann on 20/11/2021.
//

#ifndef METRICS_TYPEDEFS_H
#define METRICS_TYPEDEFS_H

#include <string>
#include <unordered_map>
#include <vector>

typedef std::vector<std::string> string_vector_t;
typedef std::vector<double> double_vector_t;
typedef std::vector<double_vector_t> double_matrix_t;
typedef std::unordered_map<char, size_t> token_index_map_t;

#endif //METRICS_TYPEDEFS_H
