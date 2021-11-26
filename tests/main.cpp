#include <iostream>
#include "PairwiseDistanceComputer.h"
#include "metrics/CdrDist.h"

int main() {
    std::string aminoAcids{"'ARNDCQEGHILKMFPSTWYVBZX*'"};

    size_t N = aminoAcids.size();
    doubleMatrix matrix (N, doubleVector (N, 0));
    stringIndexMap index;
    for (size_t i = 0; i < N; i++) {
        index[aminoAcids.substr(i, 1)] = i;
    }

    metric::CdrDist metric {matrix, index};
    stringVector inputs {std::string("CASTPGTGGLYTF"),
                         std::string("CASSHTGPLMNTEAFF"),
                         std::string("CASSQGGAINYGYTF"),
                         std::string("CASSDFGGAYNEQFF"),
                         std::string("CASSYQDRVNSPLHF")};

    PairwiseDistanceComputer executor { &metric };
    doubleVector distances = executor.computeDistance(inputs);

    std::cout << "[";
    for (const auto& distance : distances) {
        std::cout << distance;
        if (&distance == &distances.back()) break;

        std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    return 0;
}
