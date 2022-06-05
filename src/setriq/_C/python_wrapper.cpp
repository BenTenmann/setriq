//
// Created by Benjamin Tenmann on 22/11/2021.
//

#ifndef EXTENSION_NAME
#define EXTENSION_NAME _C
#endif

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pairwise_distance_computation.h"
#include "metrics/CdrDist.h"
#include "metrics/Levenshtein.h"
#include "metrics/TcrDist.h"
#include "metrics/Hamming.h"
#include "metrics/Jaro.h"
#include "metrics/JaroWinkler.h"
#include "utils/type_defs.h"

namespace py = pybind11;

// ----- pairwise distances ----------------------------------------------------------------------------------------- //
py::list cdr_dist(const string_vector_t& sequences,
                  const double_matrix_t& substitution_matrix,
                  const token_index_map_t& index,
                  const double& gap_opening_penalty,
                  const double& gap_extension_penalty) {
    metric::CdrDist metric {substitution_matrix, index, gap_opening_penalty, gap_extension_penalty};

    double_vector_t out = pairwise_distance_computation(metric, sequences);
    return py::cast(out);
}

py::list levenshtein(const string_vector_t& sequences, const double& extra_cost) {
    metric::Levenshtein metric {extra_cost};

    double_vector_t out = pairwise_distance_computation(metric, sequences);
    return py::cast(out);
}

py::list tcr_dist_component(const string_vector_t& sequences,
                            const double_matrix_t& substitution_matrix,
                            const token_index_map_t& index,
                            const double& gap_penalty,
                            const char& gap_symbol,
                            const double& distance_weight) {
    metric::TcrDist metric {substitution_matrix, index, gap_penalty, gap_symbol, distance_weight};

    double_vector_t out = pairwise_distance_computation(metric, sequences);
    return py::cast(out);
}

py::list hamming(const string_vector_t& sequences, const double& mismatch_score) {
    metric::Hamming metric {mismatch_score};

    double_vector_t out = pairwise_distance_computation(metric, sequences);
    return py::cast(out);
}

py::list jaro(const string_vector_t& sequences, const jaro_weighting_t& jaro_weights) {
    metric::Jaro metric {jaro_weights};

    double_vector_t out = pairwise_distance_computation(metric, sequences);
    return py::cast(out);
}

py::list jaro_winkler(const string_vector_t& sequences,
                      const double& p,
                      const size_t& max_l,
                      const jaro_weighting_t& jaro_weights) {
    metric::JaroWinkler metric {p, max_l, metric::Jaro{jaro_weights}};

    double_vector_t out = pairwise_distance_computation(metric, sequences);
    return py::cast(out);
}

// ----- single dispatch -------------------------------------------------------------------------------------------- //
py::float_ cdr_dist_sd(const std::string& a, std::string& b,
                       const double_matrix_t& substitution_matrix,
                       const token_index_map_t& index,
                       const double& gap_opening_penalty,
                       const double& gap_extension_penalty) {
    metric::CdrDist metric {substitution_matrix, index, gap_opening_penalty, gap_extension_penalty};
    double out = metric.forward(a, b);
    return py::cast(out);
}

py::float_ levenshtein_sd(const std::string& a, const std::string& b, const double& extra_cost) {
    metric::Levenshtein metric {extra_cost};
    double out = metric.forward(a, b);
    return py::cast(out);
}

py::float_ tcr_dist_component_sd(const std::string& a, const std::string& b,
                                 const double_matrix_t& substitution_matrix,
                                 const token_index_map_t& index,
                                 const double& gap_penalty,
                                 const char& gap_symbol,
                                 const double& distance_weight) {
    metric::TcrDist metric {substitution_matrix, index, gap_penalty, gap_symbol, distance_weight};
    double out = metric.forward(a, b);
    return py::cast(out);
}

py::float_ hamming_sd(const std::string& a, const std::string& b, const double& mismatch_score) {
    metric::Hamming metric {mismatch_score};
    double out = metric.forward(a, b);
    return py::cast(out);
}

py::float_ jaro_sd(const std::string& a, const std::string& b, const jaro_weighting_t& jaro_weights) {
    metric::Jaro metric {jaro_weights};
    double out = metric.forward(a, b);
    return py::cast(out);
}

py::float_ jaro_winkler_sd(const std::string& a, const std::string& b,
                           const double& p,
                           const size_t& max_l,
                           const jaro_weighting_t& jaro_weights) {
    metric::JaroWinkler metric {p, max_l, metric::Jaro{jaro_weights}};
    double out = metric.forward(a, b);
    return py::cast(out);
}

// ----- module def ------------------------------------------------------------------------------------------------- //
PYBIND11_MODULE(EXTENSION_NAME, m) {
    m.doc() = "Python module written in C++ for pairwise distance computation for sequences.";

    // pairwise
    m.def("cdr_dist", &cdr_dist, "Compute the pairwise CDR-dist metric for a set of CDR3 sequences.",
          py::arg("sequences"), py::arg("substitution_matrix"), py::arg("index"),
          py::arg("gap_opening_penalty"), py::arg("gap_extension_penalty"));

    m.def("levenshtein", &levenshtein, "Compute the pairwise Levenshtein distances for a set of sequences.",
          py::arg("sequences"), py::arg("extra_cost"));

    m.def("tcr_dist_component", &tcr_dist_component, "Compute pairwise TCR-dist for a set of TCR components.",
          py::arg("sequences"), py::arg("substitution_matrix"), py::arg("index"),
          py::arg("gap_penalty"), py::arg("gap_symbol"), py::arg("weight"));

    m.def("hamming", &hamming, "Compute pairwise Hamming distance for a set of sequences.",
          py::arg("sequences"), py::arg("mismatch_score"));

    m.def("jaro", &jaro, "Compute pairwise Jaro distance for a set of sequences.",
          py::arg("sequences"), py::arg("jaro_weights"));

    m.def("jaro_winkler", &jaro_winkler, "Compute pairwise Jaro-Winkler distance for a set of sequences.",
          py::arg("sequences"), py::arg("p"), py::arg("max_l"), py::arg("jaro_weights"));

    // single dispatch
    m.def("cdr_dist_sd", &cdr_dist_sd, "Compute the CDR-dist metric between two CDR3 sequences.",
          py::arg("a"), py::arg("b"), py::arg("substitution_matrix"), py::arg("index"),
          py::arg("gap_opening_penalty"), py::arg("gap_extension_penalty"));

    m.def("levenshtein_sd", &levenshtein_sd, "Compute the Levenshtein distance between two sequences.",
          py::arg("a"), py::arg("b"), py::arg("extra_cost"));

    m.def("tcr_dist_component_sd", &tcr_dist_component_sd, "Compute TCR-dist between two TCR components.",
          py::arg("a"), py::arg("b"), py::arg("substitution_matrix"), py::arg("index"),
          py::arg("gap_penalty"), py::arg("gap_symbol"), py::arg("weight"));

    m.def("hamming_sd", &hamming_sd, "Compute the Hamming distance between two sequences.",
          py::arg("a"), py::arg("b"), py::arg("mismatch_score"));

    m.def("jaro_sd", &jaro_sd, "Compute the Jaro distance between two sequences.",
          py::arg("a"), py::arg("b"), py::arg("jaro_weights"));

    m.def("jaro_winkler_sd", &jaro_winkler_sd, "Compute the Jaro-Winkler distance between two sequences.",
          py::arg("a"), py::arg("b"), py::arg("p"), py::arg("max_l"), py::arg("jaro_weights"));

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
