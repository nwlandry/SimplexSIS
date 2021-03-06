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
# graph parameters
exponent = 4 # power law exponent
minDegree = 67
maxDegree = 1000
n = 10000
simplexSize = 3
isDegreeCorrelated = True
degreeDistType = "power-law"
meanSimplexDegree = 30
meanDegree = 10
isRandom = True

#simulation parameters
timesteps = 100
dt = 0.1
numNodesToRestart = 0.0001
gamma = 2
betaCritFraction = 1.03
alpha = 0.0

# Epidemic parameters
initialFraction = 0.01
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
beta = betaCritFraction*betaCrit

plt.figure()
start = time.time()
for x0 in initialConditions:
    averageInfection, endState = simplexContagion.microscopicSimplexSISDynamics(A, simplexList, simplexIndices, gamma, beta, alpha, x0, timesteps, dt, numNodesToRestart)
    plt.plot(np.linspace(0, (timesteps-1)*dt, timesteps), averageInfection)

end = time.time()
print('The elapsed time is ' + str(end-start) + 's')
plt.xlabel("time", fontsize=18)
plt.ylabel(r"$\langle I\rangle$", fontsize=18)
plt.show()
