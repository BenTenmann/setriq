//
// Created by Benjamin Tenmann on 22/11/2021.
//

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "PairwiseDistanceComputer.h"
#include "metrics/CdrDist.h"
#include "metrics/TcrDist.h"
#include "utils/typeDefs.h"

namespace py = pybind11;

py::list cdr_dist(const stringVector &sequences, const doubleMatrix& substitutionMatrix, const stringIndexMap& index) {
    metric::CdrDist metric {substitutionMatrix, index};
    PairwiseDistanceComputer computer { &metric };

    doubleVector out = computer.computeDistance(sequences);
    return py::cast(out);
}

py::list tcr_dist_component(const stringVector& sequences,
                            const doubleMatrix& substitutionMatrix,
                            const stringIndexMap& index,
                            const double& gapPenalty,
                            const char& gapSymbol,
                            const double& distanceWeight) {
    metric::TcrDist metric {substitutionMatrix, index, gapPenalty, gapSymbol, distanceWeight};
    PairwiseDistanceComputer computer { &metric };

    doubleVector out = computer.computeDistance(sequences);
    return py::cast(out);
}

PYBIND11_MODULE(_C, m) {
    m.doc() = "Python module written in C++ for pairwise distance computation for sequences.";

    m.def("cdr_dist", &cdr_dist, "Compute the pairwise CDR-dist metric for a set of CDR3 sequences.",
          py::arg("sequences"), py::arg("substitution_matrix"), py::arg("index"));

    m.def("tcr_dist_component", &tcr_dist_component, "Compute pairwise TCR-dist for a set of TCR components.",
          py::arg("sequences"), py::arg("substitution_matrix"), py::arg("index"),
          py::arg("gap_penalty"), py::arg("gap_symbol"), py::arg("weight"));

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
