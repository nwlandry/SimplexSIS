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
k0 = 13.33557834
r = 4 # power law exponent
minDeg = 10
maxDeg = 30
n = 10000
simplexSize = 3
isIndependentUniform = True
degreeDistType = "power-law"
meanSimplexDegree = 20
meanDegree = 20

# Epidemic parameters
initialFraction = 0.01
x0 = np.random.choice([0, 1], size=n, p=[1-initialFraction, initialFraction])

#simulation parameters
timesteps = 6000
dt = 0.1
# length over which to average
avgLength = int(0.1*timesteps)
numBetaPts = 31
numAlphaPts = 9
startAlphaCritFraction = 0
endAlphaCritFraction = 2
startBetaCritFraction = 0
endBetaCritFraction = 1.5
numProcesses = numAlphaPts


# generate degree sequence and adjacency matrix
if degreeDistType == "uniform":
    degreeSequence = simplexUtilities.generateUniformDegreeSequence(n, minDeg, maxDeg)
elif degreeDistType == "power-law":
    degreeSequence = simplexUtilities.generatePowerLawDegreeSequence(n, k0, n, r)
elif degreeDistType == "poisson":
    degreeSequence = simplexUtilities.generatePoissonDegreeSequence(n, meanDegree)

A = simplexUtilities.generateConfigModelAdjacency(degreeSequence)

# Calculate values needed in critical value calculation
meanDegree = simplexUtilities.meanPowerOfDegree(A.sum(axis=0), 1)
meanSquaredDegree =  simplexUtilities.meanPowerOfDegree(A.sum(axis=0), 2)
meanCubedDegree =  simplexUtilities.meanPowerOfDegree(A.sum(axis=0), 3)

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
    [simplexList, simplexIndices] = simplexUtilities.generateConfigModelSimplexList(degreeSequence, simplexSize)
    # epidemic parameters
    gamma = 2
    betaCrit = meanDegree/meanSquaredDegree*gamma
    alphaCrit = (meanDegree**2)*meanCubedDegree/(meanSquaredDegree**3)*gamma

print("beta critical is " + str(betaCrit))
print("alpha critical is " + str(alphaCrit))

beta = np.concatenate([np.linspace(startBetaCritFraction*betaCrit, endBetaCritFraction*betaCrit, numBetaPts),
                      np.linspace(betaCrit*(endBetaCritFraction-(endBetaCritFraction-startBetaCritFraction)/(numBetaPts-1)), startBetaCritFraction*betaCrit, numBetaPts-1)])
alpha = np.linspace(startAlphaCritFraction*alphaCrit, endAlphaCritFraction*alphaCrit, numAlphaPts)

start = time.time()
equilibria = simplexContagion.generateSISEquilibria(A, simplexList, simplexIndices, gamma, beta, alpha, x0, timesteps, dt, avgLength)
end = time.time()
print('The elapsed time is ' + str(end-start) + 's')

with open('equilibriaData' + datetime.now().strftime("%m%d%Y-%H%M%S"), 'wb') as file:
    pickle.dump([gamma, beta, alpha, equilibria, betaCrit, alphaCrit, meanDegree, meanSquaredDegree, meanCubedDegree, meanSimplexDegree], file)
