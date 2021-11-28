#include <fstream>
#include <iostream>
#include <json/json.h>
#include "PairwiseDistanceComputer.h"
#include "metrics/CdrDist.h"

stringIndexMap getIndexMap (const Json::Value &jsonObject) {
    stringIndexMap index;
    for (auto const& key : jsonObject["index"].getMemberNames()) {
        index[key] = jsonObject["index"][key].asUInt();
    }
    return index;
}

doubleMatrix getSubstitutionMatrix (const Json::Value &jsonObject) {
    size_t N = jsonObject["substitution_matrix"].size();
    size_t M = jsonObject["substitution_matrix"][0].size();

    doubleMatrix matrix (N);
    for (int i = 0; i < N; ++i) {

        doubleVector row (M);
        for (int j = 0; j < M; ++j) {
            row[j] = jsonObject["substitution_matrix"][i][j].asFloat();
        }
        matrix[i] = row;
    }
    return matrix;
}

int main() {
    Json::Value root;
    Json::Reader reader;
    std::ifstream file ("data/blosum-62.json");

    reader.parse(file, root);

    stringIndexMap index {getIndexMap(root)};
    doubleMatrix matrix {getSubstitutionMatrix(root)};

    metric::CdrDist metric {matrix, index};
    stringVector inputs {std::string("AASQ"),
                         std::string("PASQ")};

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
