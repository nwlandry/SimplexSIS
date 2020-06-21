#include <iostream>
#include <stdlib.h>
//#include "SimplexContagion.hpp"
#include "SimplexUtilities.hpp"
#include <vector>
#include <tuple>
#include "SimplexContagion.hpp"
#include <chrono>
#include <iomanip>
using namespace std::chrono;

int main()
{
    // Network parameters
    int n = 1000;
    std::string degreeDistType = "power-law";
    double minDegree = 20;
    double maxDegree = 1000;
    double exponent = 4;
    bool isDegreeCorrelated = true;
    int simplexSize = 3;
    bool isMajorityVote = true;
    double meanSimplexDegree = 100;

    //Simulation parameters
    double dt = 0.1;
    int numTimesteps = 1000;
    double gamma = 2;
    double beta2 = 0.07;
    double beta3 = 0.1;
    double nodeFractionToRestart = 0.001;
    double fractionInfected = 0.001;

    // network variables
    std::vector<int> degreeSequence;
    std::vector<std::vector<int>> edgeList;
    std::vector<std::vector<int>> simplexList;
    std::vector<std::vector<int>> simplexIndices;

    //simulation variables
    std::vector<int> x0 = SimplexContagion::initializeInfectionVector(n, fractionInfected);
    std::vector<double> infected;

    degreeSequence = SimplexUtilities::generateDegreeSequence(n, degreeDistType, minDegree, maxDegree, exponent);

    SimplexUtilities::generateConfigurationModelEdgeList(degreeSequence, edgeList);

    /*for (int index = 0; index < edgeList.size(); index++)
    {
        std::cout << "(" << edgeList[index][0] << ", " << edgeList[index][1] << ") ";
    }*/
    if (isDegreeCorrelated)
    {
        SimplexUtilities::generateConfigurationSimplexModel(degreeSequence, simplexSize, simplexList, simplexIndices);
    }
    else
    {
        SimplexUtilities::generateUncorrelatedSimplices(n, meanSimplexDegree, simplexSize, simplexList, simplexIndices);
    }

    /*for (int index = 0; index < simplexList.size(); index++)
    {
        std::cout << "(" << simplexList[index][0] << ", " << simplexList[index][1] << ", " << simplexList[index][2] << ") ";
    }*/
    auto start = high_resolution_clock::now();
    infected = SimplexContagion::runMicroscopicSimplexSISDynamics(n, edgeList, simplexList, simplexIndices, gamma, beta2, beta3, isMajorityVote, x0, numTimesteps, dt, nodeFractionToRestart);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start)/1e6;

    // To get the value of duration use the count()
    // member function on the duration object
    std::cout << std::setprecision(3) << duration.count() << "s" << "\n";
    for (std::vector<double>::iterator it = infected.begin(); it != infected.end(); it++)
    {
        std::cout << *it << " ";
    }
    return 0;
}
