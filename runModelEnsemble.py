import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs
import random
import simplexUtilities
import simplexContagion
import pickle
from datetime import datetime
import time
import multiprocessing as mp

# graph parameters
r = 4 # power law exponent
minDegree = 6.695
maxDegree = 10
n = 10
simplexSize = 3
isIndependentUniform = True
degreeDistType = "power-law"
meanSimplexDegree = 100
meanDegree = 100

# Epidemic parameters
initialFraction = 0.01

#simulation parameters
numProcesses = 24
numSimulations = 10
timesteps = 1000
dt = 0.1
gamma = 2
# length over which to average
avgLength = int(0.3*timesteps)
numBetaPts = 41
numAlphaPts = 24
startAlpha = 0
endAlpha = 0.075
startBetaCritFraction = 0
endBetaCritFraction = 1.5

# generate degree sequence and adjacency matrix
degreeSequenceList = list()
adjacencyList = list()
for i in range(numSimulations):
    if degreeDistType == "uniform":
        degreeSequence = simplexUtilities.generateUniformDegreeSequence(n, minDegree, maxDegree)
    elif degreeDistType == "power-law":
        degreeSequence = simplexUtilities.generatePowerLawDegreeSequence(n, minDegree, maxDegree, r)
    elif degreeDistType == "poisson":
        degreeSequence = simplexUtilities.generatePoissonDegreeSequence(n, meanDegree)
    degreeSequenceList.append(degreeSequence)
    adjacencyList.append(simplexUtilities.generateConfigModelAdjacency(degreeSequence))

# Calculate values needed in critical value calculation
meanDegree = 0
meanSquaredDegree = 0
meanCubedDegree = 0
meanNumSelfLoops = 0
for i in range(numSimulations):
    meanDegree = meanDegree + simplexUtilities.meanPowerOfDegree(degreeSequenceList[i], 1)/numSimulations
    meanSquaredDegree = meanSquaredDegree + simplexUtilities.meanPowerOfDegree(degreeSequenceList[i], 2)/numSimulations
    meanCubedDegree = meanCubedDegree + simplexUtilities.meanPowerOfDegree(degreeSequenceList[i], 3)/numSimulations
    meanNumSelfLoops = meanNumSelfLoops + np.trace(adjacencyList[i].todense())/numSimulations

betaCrit = meanDegree/meanSquaredDegree*gamma
if isIndependentUniform:
    alphaCrit = meanCubedDegree/(meanDegree**3 * meanSimplexDegree)*gamma
else:
    alphaCrit = (meanDegree**2)*meanCubedDegree/(meanSquaredDegree**3)*gamma

print("The mean degree is {:.2f}".format(meanDegree))
print("The mean squared degree is {:.2f}".format(meanSquaredDegree))
print("The mean cubed degree is {:.2f}".format(meanCubedDegree))
print("The mean simplex degree is {}".format(meanSimplexDegree))
print("{} self-loops".format(meanNumSelfLoops))
print("beta critical is " + str(betaCrit))
print("alpha critical is " + str(alphaCrit))

#Generate simplex list
simplexSetList = list()
simplexIndicesList = list()

for i in range(numSimulations):
    if isIndependentUniform:
        [simplexSet, simplexIndices] = simplexUtilities.generateUniformSimplexList(n, meanSimplexDegree, simplexSize)
        simplexSetList.append(simplexSet)
        simplexIndicesList.append(simplexIndices)
    else:
        [simplexSet, simplexIndices] = simplexUtilities.generateConfigModelSimplexList(degreeSequence, simplexSize)
        simplexSetList.append(simplexSet)
        simplexIndicesList.append(simplexIndices)



beta = np.concatenate([np.linspace(startBetaCritFraction*betaCrit, endBetaCritFraction*betaCrit, numBetaPts),
                      np.linspace(betaCrit*(endBetaCritFraction-(endBetaCritFraction-startBetaCritFraction)/(numBetaPts-1)), startBetaCritFraction*betaCrit, numBetaPts-1)])
alpha = np.linspace(startAlpha, endAlpha, numAlphaPts)

start = time.time()
averagedEquilibria = simplexContagion.generateSISEquilibriaEnsembleParallelized(adjacencyList, simplexSetList, simplexIndicesList, gamma, beta, alpha, initialFraction, timesteps, dt, avgLength, numSimulations, numProcesses)
end = time.time()
print('The elapsed time is ' + str(end-start) + 's')

with open('equilibriaData' + datetime.now().strftime("%m%d%Y-%H%M%S"), 'wb') as file:
    pickle.dump([gamma, beta, alpha, averagedEquilibria, equilibria, betaCrit, meanDegree, meanSquaredDegree, meanCubedDegree, meanSimplexDegree], file)