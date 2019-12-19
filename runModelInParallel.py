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
k0 = 50
r = 4 # power law exponent
minDeg = 50
maxDeg = 100
n = 10000
simplexSize = 3
isIndependentUniform = True
degreeDistType = "poisson"

# Epidemic parameters
initialFraction = 0.5
x0 = np.random.choice([0, 1], size=n, p=[1-initialFraction, initialFraction])

#simulation parameters
timesteps = 6000
dt = 1
# length over which to average
avgLength = int(0.1*timesteps)
numBetaPts = 31
numAlphaPts = 8
startAlphaCritFraction = 0
endAlphaCritFraction = 4
startBetaCritFraction = 0.5
endBetaCritFraction = 1.5
numProcesses = numAlphaPts
meanDegree = 20


# generate degree sequence and adjacency matrix
if degreeDistType == "uniform":
    k = simplexUtilities.generateUniformDegreeSequence(n, minDeg, maxDeg)
elif degreeDistType == "power-law":
    k = simplexUtilities.generatePowerLawDegreeSequence(n, k0, n, r)
elif degreeDistType == "poisson":
    k = simplexUtilities.generatePoissonDegreeSequence(n, meanDegree)

A = simplexUtilities.generateConfigModelAdjacency(k)

# Calculate critical values
meanDegree = simplexUtilities.meanPowerOfDegree(A.sum(axis=0), 1)
meanSquaredDegree =  simplexUtilities.meanPowerOfDegree(A.sum(axis=0), 2)
meanCubedDegree =  simplexUtilities.meanPowerOfDegree(A.sum(axis=0), 3)
meanSimplexDegree = 6

print("The mean degree is {:.2f}".format(meanDegree))
print("The mean squared degree is {:.2f}".format(meanSquaredDegree))
print("The mean cubed degree is {:.2f}".format(meanCubedDegree))
print("The mean simplex degree is {}".format(meanSimplexDegree))

print("{} self-loops".format(np.trace(A.todense())))


#Generate simplex list
if isIndependentUniform:
    #[simplexList, simplexIndices] = simplexUtilities.generateUniformSimplexList(n, int(meanDegree*n), simplexSize)
    [simplexList, simplexIndices] = simplexUtilities.generateUniformSimplexList(n, meanSimplexDegree, simplexSize)
    # epidemic parameters
    gamma = 2
    betaCrit = meanDegree/meanSquaredDegree*gamma
    alphaCrit = meanCubedDegree/(meanDegree**3 * meanSimplexDegree)*gamma

else:
    [simplexList, simplexIndices] = simplexUtilities.generateConfigModelSimplexList(k, simplexSize)
    # epidemic parameters
    gamma = 2
    betaCrit = meanDegree/meanSquaredDegree*gamma
    alphaCrit = (meanDegree**2)*meanCubedDegree/(meanSquaredDegree**3)*gamma

print("beta critical is " + str(betaCrit))
print("alpha critical is " + str(alphaCrit))

# beta = np.concatenate([np.linspace(startBetaCritFraction*betaCrit, endBetaCritFraction*betaCrit, numBetaPts),
#                       np.linspace(betaCrit*(endBetaCritFraction-(endBetaCritFraction-startBetaCritFraction)/(numBetaPts-1)), startBetaCritFraction*betaCrit, numBetaPts-1)])
# alpha = np.linspace(startAlphaCritFraction*alphaCrit, endAlphaCritFraction*alphaCrit, numAlphaPts)
#
# start = time.time()
# equilibria = simplexContagion.generateSISEquilibriaParallelized(A, simplexList, simplexIndices, gamma, beta, alpha, x0, timesteps, dt, avgLength, numProcesses)
# end = time.time()
# print('The elapsed time is ' + str(end-start) + 's')
#
# with open('equilibriaData' + datetime.now().strftime("%m%d%Y-%H%M%S"), 'wb') as file:
#     pickle.dump([alpha, beta, gamma, meanDegree, meanSimplexDegree, equilibria], file)
