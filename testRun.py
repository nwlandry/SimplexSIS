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
import simplexTheory

# graph parameters
# graph parameters
exponent = 3 # power law exponent
minDegree = 10
maxDegree = 30
n = 30000
simplexSize = 3
isDegreeCorrelated = True
degreeDistType = "uniform"
meanSimplexDegree = 20
meanDegree = 20
isRandom = True

#simulation parameters
timesteps = 3000
dt = 0.1
numNodesToRestart = 10
gamma = 2
betaCritFraction = 0.99
alphaCritFraction = 1.5
tolerance = 0.0001
minAlpha = 0
maxAlpha = 0.5
digits = 5
option = "fast"

# Epidemic parameters
initialFraction = 0.0
x01 = np.random.choice([0, 1], size=n, p=[1-initialFraction, initialFraction])
# initialFraction = 1
# x02 = np.random.choice([0, 1], size=n, p=[1-initialFraction, initialFraction])
initialConditions = [x01]

# generate degree sequence and adjacency matrix
if degreeDistType == "uniform":
    degreeSequence = simplexUtilities.generateUniformDegreeSequence(n, minDegree, maxDegree, isRandom=isRandom)
elif degreeDistType == "power-law":
    degreeSequence = simplexUtilities.generatePowerLawDegreeSequence(n, minDegree, maxDegree, exponent, isRandom=isRandom)
elif degreeDistType == "poisson":
    degreeSequence = simplexUtilities.generatePoissonDegreeSequence(n, meanDegree)
degreeHist = simplexTheory.degreeSequenceToHist(degreeSequence)

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
alphaCrit = simplexTheory.calculateTheoreticalCriticalAlpha(gamma, betaCrit, minAlpha, maxAlpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits, tolerance=tolerance, option=option)
print(stop)
alpha = alphaCritFraction*alphaCrit
beta = betaCritFraction*betaCrit

roots = simplexTheory.solveEquilibrium(gamma, beta, alpha, degreeHist, meanSimplexDegree=None, isDegreeCorrelated=isDegreeCorrelated, majorityVote=True, healing=False, digits=4)

plt.figure()
for root in roots:
    plt.plot([0, (timesteps-1)*dt], [root, root], 'k--')

start = time.time()
for x0 in initialConditions:
    averageInfection, endState = simplexContagion.microscopicSimplexSISDynamics(A, simplexList, simplexIndices, gamma, beta, alpha, x0, timesteps, dt, numNodesToRestart)
    plt.plot(np.linspace(0, (timesteps-1)*dt, timesteps), averageInfection)

end = time.time()
print('The elapsed time is ' + str(end-start) + 's')
plt.xlabel("time", fontsize=18)
plt.ylabel(r"$\langle I\rangle$", fontsize=18)
plt.show()
