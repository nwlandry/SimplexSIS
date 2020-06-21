#pragma once
#include <vector>

namespace SimplexContagion
{
    std::vector<double> runMicroscopicSimplexSISDynamics(int n, std::vector<std::vector<int>> edgeList, std::vector<std::vector<int>> simplexList, std::vector<std::vector<int>> simplexIndices, double gamma, double beta2, double beta3, bool isMajorityVote, std::vector<int> x0, int numTimesteps, double dt, double nodeFractionToRestart);

    std::vector<double> calculateNumberOfInfectedNeighbors(std::vector<std::vector<int>> edgeList, std::vector<int> x);

    std::vector<double> calculateInfectionProbabilitiesByNeighbors(std::vector<double> numberOfInfectedNeighbors, double beta2, double dt);

    double calculateNumberOfInfectedSimplexNeighbors(int index, std::vector<std::vector<int>> simplexIndices, std::vector<std::vector<int>> simplexList, std::vector<int> x, bool isMajorityVote);

    double calculateInfectionProbabilityBySimplexNeighbors(double numberOfInfectedSimplexNeighbors, double beta3, double dt);

    std::vector<int> initializeInfectionVector(int n, double fractionInfected);
};
