//
// Created by Benjamin Tenmann on 22/11/2021.
//

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "PairwiseDistanceComputer.h"
#include "metrics/CdrDist.h"
#include "utils/typeDefs.h"

namespace py = pybind11;

py::list cdr_dist(const stringVector &sequences, const doubleMatrix& substitutionMatrix, const stringIndexMap& index) {
    metric::CdrDist metric {substitutionMatrix, index};
    PairwiseDistanceComputer computer { &metric };

    doubleVector out = computer.computeDistance(sequences);
    return py::cast(out);
}

PYBIND11_MODULE(setriq, m) {
    m.doc() = "Python module written in C++ for pairwise distance computation for sequences.";

    m.def("cdr_dist", &cdr_dist, "Compute the pairwise CDR-dist metric for a set of CDR3 sequences.",
          py::arg("sequences"), py::arg("substitution_matrix"), py::arg("index"));

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
