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
minDegree = 66.68
maxDegree = 10000
n = 10000
simplexSize = 3
isIndependent = True
degreeDistType = "power-law"
meanSimplexDegree = 100
meanDegree = 100
isRandom = False

# Epidemic parameters
initialFraction = 0.01
x0 = np.random.choice([0, 1], size=n, p=[1-initialFraction, initialFraction])

#simulation parameters
timesteps = 1000
dt = 0.1
# length over which to average
avgLength = int(0.3*timesteps)
numBetaPts = 31
numAlphaPts = 24
gamma = 2
alphaCrit = 0.066
startAlphaCritFraction = 0.5
endAlphaCritFraction = 1.5
startBetaCritFraction = 0.5
endBetaCritFraction = 1.5
numProcesses = numAlphaPts


# generate degree sequence and adjacency matrix
if degreeDistType == "uniform":
    degreeSequence = simplexUtilities.generateUniformDegreeSequence(n, minDegree, maxDegree, isRandom=isRandom)
elif degreeDistType == "power-law":
    degreeSequence = simplexUtilities.generatePowerLawDegreeSequence(n, minDegree, maxDegree, r, isRandom=isRandom)
elif degreeDistType == "poisson":
    degreeSequence = simplexUtilities.generatePoissonDegreeSequence(n, meanDegree)

A = simplexUtilities.generateConfigModelAdjacency(degreeSequence)

# Calculate values needed in critical value calculation
meanDegree = simplexUtilities.meanPowerOfDegree(degreeSequence, 1)
meanSquaredDegree =  simplexUtilities.meanPowerOfDegree(degreeSequence, 2)
meanCubedDegree =  simplexUtilities.meanPowerOfDegree(degreeSequence, 3)

print("The mean degree is {:.2f}".format(meanDegree))
print("The mean squared degree is {:.2f}".format(meanSquaredDegree))
print("The mean cubed degree is {:.2f}".format(meanCubedDegree))
print("{} self-loops".format(np.trace(A.todense())))


#Generate simplex list
if isIndependent:
    [simplexList, simplexIndices] = simplexUtilities.generateUniformSimplexList(n, meanSimplexDegree, simplexSize)
else:
    [simplexList, simplexIndices] = simplexUtilities.generateConfigModelSimplexList(degreeSequence, simplexSize)

betaCrit = meanDegree/meanSquaredDegree*gamma
meanSimplexDegree = len(simplexList)/n
print("The mean simplex degree is {}".format(meanSimplexDegree))

print("beta critical is " + str(betaCrit))

beta = np.concatenate([np.linspace(startBetaCritFraction*betaCrit, endBetaCritFraction*betaCrit, numBetaPts),
                      np.linspace(betaCrit*(endBetaCritFraction-(endBetaCritFraction-startBetaCritFraction)/(numBetaPts-1)), startBetaCritFraction*betaCrit, numBetaPts-1)])
alpha = np.linspace(startAlphaCritFraction*alphaCrit, endAlphaCritFraction*alphaCrit, numAlphaPts)

start = time.time()
equilibria = simplexContagion.generateSISEquilibriaParallelized(A, simplexList, simplexIndices, gamma, beta, alpha, x0, timesteps, dt, avgLength, numProcesses)
end = time.time()
print('The elapsed time is ' + str(end-start) + 's')

with open('equilibriaData' + datetime.now().strftime("%m%d%Y-%H%M%S"), 'wb') as file:
    pickle.dump([gamma, beta, alpha, equilibria, degreeSequence, isIndependent, degreeDistType, r, meanSimplexDegree], file)
