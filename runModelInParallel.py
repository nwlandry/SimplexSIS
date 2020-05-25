import random
import simplexUtilities
import simplexContagion
import pickle
from datetime import datetime
import time
import simplexTheory
import os
import sys
import numpy as np

# graph parameters
degreeDistType = sys.argv[1]
minDegree = float(sys.argv[2])
maxDegree = float(sys.argv[3])
n = int(sys.argv[4])
isDegreeCorrelated = (sys.argv[5] == "True")
meanSimplexDegree = float(sys.argv[6])
exponent = float(sys.argv[7]) # power law exponent
isRandom = True
simplexSize = 3
meanDegree = 100

# Epidemic parameters
initialFraction = 0.01
x0 = np.random.choice([0, 1], size=n, p=[1-initialFraction, initialFraction])

#simulation parameters
numProcesses = len(os.sched_getaffinity(0))
print("Number of cores is " + str(numProcesses))
timesteps = 1000
dt = 0.1
nodeFractionToRestart = 0.0002
# length over which to average
avgLength = int(0.5*timesteps)
numBetaPts = 31
numAlphaPts = numProcesses
gamma = 2

tolerance = 0.0001
option = "fast"
minAlpha = 0
maxAlpha = 0.2
digits = 5

startAlphaCritFraction = 0.5
endAlphaCritFraction = 1.5
startBetaCritFraction = 0.5
endBetaCritFraction = 1.5

# generate degree sequence and adjacency matrix
if degreeDistType == "uniform":
    degreeSequence = simplexUtilities.generateUniformDegreeSequence(n, minDegree, maxDegree, isRandom=isRandom)
elif degreeDistType == "power-law":
    degreeSequence = simplexUtilities.generatePowerLawDegreeSequence(n, minDegree, maxDegree, exponent, isRandom=isRandom)
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
if isDegreeCorrelated:
    [simplexList, simplexIndices] = simplexUtilities.generateConfigModelSimplexList(degreeSequence, simplexSize)
else:
    [simplexList, simplexIndices] = simplexUtilities.generateUniformSimplexList(n, meanSimplexDegree, simplexSize)

betaCrit = meanDegree/meanSquaredDegree*gamma
meanSimplexDegree = simplexSize*len(simplexList)/n

degreeHist = simplexTheory.degreeSequenceToHist(degreeSequence)

alphaCrit = simplexTheory.calculateTheoreticalCriticalAlpha(gamma, betaCrit, minAlpha, maxAlpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits, tolerance=tolerance, option=option)

print("The mean simplex degree is {}".format(meanSimplexDegree))
print("beta critical is " + str(betaCrit))
print("alpha critical is " + str(alphaCrit))

beta = np.concatenate([np.linspace(startBetaCritFraction*betaCrit, endBetaCritFraction*betaCrit, numBetaPts),
                      np.linspace(betaCrit*(endBetaCritFraction-(endBetaCritFraction-startBetaCritFraction)/(numBetaPts-1)), startBetaCritFraction*betaCrit, numBetaPts-1)])
alpha = np.linspace(startAlphaCritFraction*alphaCrit, endAlphaCritFraction*alphaCrit, numAlphaPts)

start = time.time()
equilibria = simplexContagion.generateSISEquilibriaParallelized(A, simplexList, simplexIndices, gamma, beta, alpha, x0, timesteps, dt, avgLength, nodeFractionToRestart, numProcesses)
end = time.time()
print('The elapsed time is ' + str(end-start) + 's')

with open('equilibriaData' + datetime.now().strftime("%m%d%Y-%H%M%S"), 'wb') as file:
    pickle.dump([gamma, beta, alpha, equilibria, degreeSequence, isDegreeCorrelated, degreeDistType, exponent, meanSimplexDegree], file)
