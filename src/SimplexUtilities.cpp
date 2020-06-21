#include "SimplexUtilities.hpp"
#include <algorithm>
#include <iterator>
#include <random>
#include <iostream>
#include <stdlib.h>
#include <math.h>

std::vector<int> SimplexUtilities::generateDegreeSequence(int n, std::string degreeDistType, double minDegree, double maxDegree, double exponent)
{
    std::vector<int> degreeSequence;
    for (int step = 0; step < n; step++)
    {
        degreeSequence.push_back(computeRandomDegree(degreeDistType, minDegree, maxDegree, exponent));
    }
    //std::sort(degreeSequence.begin(), degreeSequence.end());
    return degreeSequence;
}

int SimplexUtilities::computeRandomDegree(std::string degreeDistType, double minDegree, double maxDegree, double exponent)
{
    double u = std::rand()/(double)RAND_MAX;
    if (degreeDistType.compare("power-law"))
    {
        return (int) minDegree*(1 - u) + maxDegree*u;
    }
    else if (degreeDistType.compare("uniform"))
    {
        return (int) powerLawDegree(u, minDegree, maxDegree, exponent);
    }
    else
    {
        return 0;
    }
}

int SimplexUtilities::powerLawDegree(double u, double minDegree, double maxDegree, double exponent)
{
    return pow(pow(minDegree, (1-exponent)) + u*(pow(maxDegree, (1-exponent)) - pow(minDegree, (1-exponent))), 1/(1-exponent));
}

void SimplexUtilities::generateConfigurationModelEdgeList(std::vector<int> degreeSequence, std::vector<std::vector<int>> &edgeListInput)
{
    std::vector<int> stubs;
    std::vector<std::vector<int>> edgeList;

    std::random_device randDevice;
    unsigned seed = randDevice();
    // Initialize a default_random_engine with the seed
    std::default_random_engine randEngine(seed);

    // If an odd number of degrees
    int n = degreeSequence.size();
    int remainder = std::accumulate(degreeSequence.begin(), degreeSequence.end(), 0) % 2;
    if (remainder != 0)
    {
        int i = std::rand() % n;
        degreeSequence[i]++;
    }
    // Generate stubs
    for (int index = 0; index < n; index++)
    {
        for (int numStub = 0; numStub < degreeSequence[index]; numStub++)
        {
            stubs.push_back(index);
        }
    }
    // Shuffle only once
    std::shuffle(stubs.begin(), stubs.end(), randEngine);
    //Assign edges
    while (stubs.size() != 0)
    {
        // Shuffle the stubs
        //std::shuffle(stubs.begin(), stubs.end(), randEngine);
        // Get the first two values
        int i = stubs[0];
        int j = stubs[1];
        edgeList.push_back(std::vector<int>{i, j});
        stubs.erase(stubs.begin(), stubs.begin() + 2);
    }
    edgeListInput = edgeList;
}


void SimplexUtilities::generateConfigurationSimplexModel(std::vector<int> degreeSequence, int simplexSize, std::vector<std::vector<int>> &simplexListInput, std::vector<std::vector<int>> &simplexIndicesInput)
{
    std::vector<int> stubs;
    std::random_device randDevice;
    unsigned seed = randDevice();
    int n = degreeSequence.size();
    std::vector<int> indices(simplexSize);
    std::vector<std::vector<int>> simplexIndices(n);
    std::vector<std::vector<int>> simplexList;
    // Initialize a default_random_engine with the seed
    std::default_random_engine randEngine(seed);

    // If an odd number of degrees
    int remainder = std::accumulate(degreeSequence.begin(), degreeSequence.end(), 0) % simplexSize;
    if (remainder != 0)
    {
        for (int addEdges = 0; addEdges < simplexSize - remainder; addEdges++)
        {
            int i = std::rand() % n;
            degreeSequence[i]++;
        }
    }
    // Generate stubs
    for (int index = 0; index < n; index++)
    {
        for (int numStub = 0; numStub < degreeSequence[index]; numStub++)
        {
            stubs.push_back(index);
        }
    }
    std::shuffle(stubs.begin(), stubs.end(), randEngine);
    //Assign edges
    while (stubs.size() != 0)
    {
        // Shuffle the stubs
        //std::shuffle(stubs.begin(), stubs.end(), randEngine);
        // Get the first two values
        for (int i = 0; i < simplexSize; i++)
        {
            simplexIndices[stubs[i]].push_back(simplexList.size()); // This is the index because we are 0 indexing.
            indices[i] = stubs[i];
        }
        simplexList.push_back(indices);
        stubs.erase(stubs.begin(), stubs.begin() + simplexSize);
    }
    simplexIndicesInput = simplexIndices;
    simplexListInput = simplexList;
}

void SimplexUtilities::generateUncorrelatedSimplices(int n, double meanSimplexDegree, int simplexSize, std::vector<std::vector<int>> &simplexListInput, std::vector<std::vector<int>> &simplexIndicesInput)
{
    std::random_device randDevice;
    unsigned seed = randDevice();
    int m = (int) meanSimplexDegree*n/simplexSize;
    std::vector<int> indices(simplexSize);
    std::vector<std::vector<int>> simplexIndices(n);
    std::vector<std::vector<int>> simplexList(m);
    // Initialize a default_random_engine with the seed
    std::default_random_engine randEngine(seed);
    for (int index = 0; index < m; index++)
    {
        for (int node = 0; node < simplexSize; node++)
        {
            int i = std::rand() % n;
            indices[node] = i;
            simplexIndices[i].push_back(simplexList.size());
        }
        simplexList.push_back(indices);
    }
    simplexListInput = simplexList;
    simplexIndicesInput = simplexIndices;
}
