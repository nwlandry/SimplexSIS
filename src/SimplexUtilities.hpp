#include <vector>
#include <string>

namespace SimplexUtilities
{
    std::vector<int> generateDegreeSequence(int n, std::string degreeDistType, double minDegree, double maxDegree, double exponent);

    int computeRandomDegree(std::string degreeDistType, double minDegree, double maxDegree, double exponent);

    int powerLawDegree(double u, double minDegree, double maxDegree, double exponent);

    void generateConfigurationModelEdgeList(std::vector<int> degreeSequence, std::vector<std::vector<int>> &edgeList);

    void generateConfigurationSimplexModel(std::vector<int> degreeSequence, int simplexSize, std::vector<std::vector<int>> &simplexList, std::vector<std::vector<int>> &simplexIndices);

    void generateUncorrelatedSimplices(int n, double meanSimplexDegree, int simplexSize, std::vector<std::vector<int>> &simplexListInput, std::vector<std::vector<int>> &simplexIndicesInput);
};
