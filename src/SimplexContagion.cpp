#include "SimplexContagion.hpp"
#include <math.h>
#include <random>
#include <tuple>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>

std::vector<double> SimplexContagion::runMicroscopicSimplexSISDynamics(int n, std::vector<std::vector<int>> edgeList, std::vector<std::vector<int>> simplexList, std::vector<std::vector<int>> simplexIndices, double gamma, double beta2, double beta3, bool isMajorityVote, std::vector<int> x0, int numTimesteps, double dt, double nodeFractionToRestart)
{
    std::vector<double> averageInfected;
    std::vector<int> x = x0;
    std::vector<int> xNew = x;
    std::vector<double> pInfectByNeighbors;

    for (int step = 0; step < numTimesteps; step++)
    {
        //std::cout << "step " << step << "\n" << std::flush;
        pInfectByNeighbors = calculateInfectionProbabilitiesByNeighbors(calculateNumberOfInfectedNeighbors(edgeList, x), beta2, dt);
        // restart nodes if infection dies out
        if (std::accumulate(x.begin(),x.end(), 0) == 0)
        {
            for (int restart = 0; restart < int(n*nodeFractionToRestart); restart++)
            {
                int i = std::rand() % n;
                if (x[i] == 0)
                {
                    x[i] = 1;
                }
            }
        }
        for (int index = 0 ; index < x.size(); index++)
        {
            // If healthy
            if (x[index] == 1)
            {
                //double u = std::rand()/(double)RAND_MAX;
                xNew[index] = (std::rand()/(double)RAND_MAX < gamma*dt ? 0 : 1);
            }
            // If infected
            else
            {
                //double u = std::rand()/(double)RAND_MAX;
                xNew[index] = (std::rand()/(double)RAND_MAX < pInfectByNeighbors[index] ? 1 : 0);
                if (xNew[index] == 0 && beta3 != 0)
                {
                    //double u = std::rand()/(double)RAND_MAX;
                    xNew[index] = (std::rand()/(double)RAND_MAX < calculateInfectionProbabilityBySimplexNeighbors( calculateNumberOfInfectedSimplexNeighbors(index, simplexList, simplexIndices, x, isMajorityVote), beta3, dt) ? 1 : 0);
                }
            }
        }
        averageInfected.push_back(((double) std::accumulate(x.begin(),x.end(), 0))/((double) x.size()));
        x = xNew;
    }
    //averageInfected.push_back(((double) std::accumulate(x.begin(),x.end(), 0))/((double) x.size()));
    return averageInfected;
}

std::vector<double> SimplexContagion::calculateNumberOfInfectedNeighbors(std::vector<std::vector<int>> edgeList, std::vector<int> x)
{
    int n = x.size();
    std::vector<double> numberOfInfectedNeighbors(n, 0.0);

    for (int index = 0 ; index < edgeList.size(); index++)
    {
        numberOfInfectedNeighbors[edgeList[index][0]] = numberOfInfectedNeighbors[edgeList[index][0]] + x[edgeList[index][1]];
        numberOfInfectedNeighbors[edgeList[index][1]] = numberOfInfectedNeighbors[edgeList[index][1]] + x[edgeList[index][0]];
    }
    return numberOfInfectedNeighbors;
}

std::vector<double> SimplexContagion::calculateInfectionProbabilitiesByNeighbors(std::vector<double> numberOfInfectedNeighbors, double beta2, double dt)
{
    std::vector<double> pInfectByNeighbors;

    for (std::vector<double>::iterator it = numberOfInfectedNeighbors.begin() ; it != numberOfInfectedNeighbors.end(); it++)
    {
        pInfectByNeighbors.push_back(1.0 - pow((1.0 - beta2*dt), *it));
    }
    return pInfectByNeighbors;
}

double SimplexContagion::calculateNumberOfInfectedSimplexNeighbors(int index, std::vector<std::vector<int>> simplexList, std::vector<std::vector<int>> simplexIndices, std::vector<int> x, bool isMajorityVote)
{
    int n = x.size();
    double numberOfInfectedSimplexNeighbors = 0.0;
    int threshold = (isMajorityVote ? 2 : 1);
    int sum = 0;

    for (std::vector<int>::iterator it = simplexIndices[index].begin() ; it != simplexIndices[index].end(); it++)
    {
        sum = 0;
        for (std::vector<int>::iterator index = simplexList[*it].begin(); index != simplexList[*it].end(); index++)
        {
            sum += x[*index];
        }
        if (sum >= threshold)
        {
            numberOfInfectedSimplexNeighbors = numberOfInfectedSimplexNeighbors + 1.0;
        }

    }
    return numberOfInfectedSimplexNeighbors;
}

double SimplexContagion::calculateInfectionProbabilityBySimplexNeighbors(double numberOfInfectedSimplexNeighbors, double beta3, double dt)
{
    return 1.0 - pow((1.0 - beta3*dt), numberOfInfectedSimplexNeighbors);
}

std::vector<int> SimplexContagion::initializeInfectionVector(int n, double fractionInfected)
{
    std::vector<int> infected(n, 0);
    for (std::vector<int>::iterator it = infected.begin(); it != infected.end(); it++)
    {
        *it = (std::rand()/(double)RAND_MAX < fractionInfected ? 1 : 0);
    }
    return infected;
}
