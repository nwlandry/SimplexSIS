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
isIndependentUniform = False
degreeDistType = "power-law"

# Epidemic parameters
initialFraction = 0.5
x0 = np.random.choice([0, 1], size=n, p=[1-initialFraction, initialFraction])

#simulation parameters
timesteps = 1000
dt = 0.1
avgLength = int(0.4*timesteps)
numBetaPts = 31
numAlphaPts = 8
startAlphaCritFraction = 0.5
endAlphaCritFraction = 4
startBetaCritFraction = 0.5
endBetaCritFraction = 1.5
numProcesses = numAlphaPts


# generate degree sequence and adjacency matrix
if degreeDistType == "uniform":
    k = simplexUtilities.generateUniformDegreeSequence(n, minDeg, maxDeg)
elif degreeDistType == "power-law":
    k = simplexUtilities.generatePowerLawDegreeSequence(n, k0, n, r)

A = simplexUtilities.generateConfigModelAdjacency(k)

# Calculate critical values
kAvg = simplexUtilities.avgPowerK(A.sum(axis=0), 1)
kSquaredAvg =  simplexUtilities.avgPowerK(A.sum(axis=0), 2)
kCubedAvg =  simplexUtilities.avgPowerK(A.sum(axis=0), 3)
print(kAvg)
print(kSquaredAvg)
print(kCubedAvg)

#Generate simplex list
if isIndependentUniform:
    [simplexList, simplexIndices] = simplexUtilities.generateUniformSimplexList(n, int(kAvg)*n, simplexSize)
    # epidemic parameters
    gamma = 2
    betaCrit = kAvg/kSquaredAvg*gamma
    alphaCrit = kCubedAvg/(kAvg**4)*gamma

else:
    [simplexList, simplexIndices] = simplexUtilities.generateConfigModelSimplexList(k, simplexSize)
    # epidemic parameters
    gamma = 2
    betaCrit = kAvg/kSquaredAvg*gamma
    alphaCrit = (kAvg**2)*kCubedAvg/(kSquaredAvg**3)*gamma



print("beta critical is " + str(betaCrit))
print("alpha critical is " + str(alphaCrit))

beta = np.concatenate([np.linspace(startBetaCritFraction*betaCrit, endBetaCritFraction*betaCrit, numBetaPts),
                      np.linspace(betaCrit*(endBetaCritFraction-(endBetaCritFraction-startBetaCritFraction)/(numBetaPts-1)), startBetaCritFraction*betaCrit, numBetaPts-1)])
alpha = np.linspace(startAlphaCritFraction*alphaCrit, endAlphaCritFraction*alphaCrit, numAlphaPts)

start = time.time()
equilibria = simplexContagion.generateSISEquilibriaParallelized(A, simplexList, simplexIndices, gamma, beta, alpha, x0, timesteps, dt, avgLength, numProcesses)
end = time.time()
print('The elapsed time is ' + str(end-start) + 's')

with open('equilibriaData' + datetime.now().strftime("%m%d%Y-%H%M%S"), 'wb') as file:
    pickle.dump([alpha, beta, equilibria], file)
