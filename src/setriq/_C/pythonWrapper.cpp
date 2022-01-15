//
// Created by Benjamin Tenmann on 22/11/2021.
//

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "PairwiseDistanceComputer.h"
#include "metrics/CdrDist.h"
#include "metrics/Levenshtein.h"
#include "metrics/TcrDist.h"
#include "utils/typeDefs.h"

namespace py = pybind11;

py::list cdr_dist(const string_vector_t& sequences,
                  const double_matrix_t& substitution_matrix,
                  const token_index_map_t& index,
                  const double& gap_opening_penalty,
                  const double& gap_extension_penalty) {
    metric::CdrDist metric {substitution_matrix, index, gap_opening_penalty, gap_extension_penalty};
    PairwiseDistanceComputer computer { &metric };

    double_vector_t out = computer.compute_distance(sequences);
    return py::cast(out);
}

py::list levenshtein(const string_vector_t& sequences, const double& extra_cost) {
    metric::Levenshtein metric {extra_cost};
    PairwiseDistanceComputer computer { &metric };

    double_vector_t out = computer.compute_distance(sequences);
    return py::cast(out);
}

py::list tcr_dist_component(const string_vector_t& sequences,
                            const double_matrix_t& substitution_matrix,
                            const token_index_map_t& index,
                            const double& gap_penalty,
                            const char& gap_symbol,
                            const double& distance_weight) {
    metric::TcrDist metric {substitution_matrix, index, gap_penalty, gap_symbol, distance_weight};
    PairwiseDistanceComputer computer { &metric };

    double_vector_t out = computer.compute_distance(sequences);
    return py::cast(out);
}

PYBIND11_MODULE(_C, m) {
    m.doc() = "Python module written in C++ for pairwise distance computation for sequences.";

    m.def("cdr_dist", &cdr_dist, "Compute the pairwise CDR-dist metric for a set of CDR3 sequences.",
          py::arg("sequences"), py::arg("substitution_matrix"), py::arg("index"),
          py::arg("gap_opening_penalty"), py::arg("gap_extension_penalty"));

    m.def("levenshtein", &levenshtein, "Compute the pairwise Levenshtein distances for a set of sequences.",
          py::arg("sequences"), py::arg("extra_cost"));

    m.def("tcr_dist_component", &tcr_dist_component, "Compute pairwise TCR-dist for a set of TCR components.",
          py::arg("sequences"), py::arg("substitution_matrix"), py::arg("index"),
          py::arg("gap_penalty"), py::arg("gap_symbol"), py::arg("weight"));

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
