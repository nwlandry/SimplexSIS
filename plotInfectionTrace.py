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
maxDegree = 1000
n = 1000
simplexSize = 3
isIndependentUniform = True
degreeDistType = "power-law"
meanSimplexDegree = 100
meanDegree = 100


# Epidemic parameters
initialFraction = 0.01
x0 = np.random.choice([0, 1], size=n, p=[1-initialFraction, initialFraction])

#simulation parameters
timesteps = 1000
dt = 0.1
alpha = 0.04
betaCritFraction = 1.1
gamma = 2

# generate degree sequence and adjacency matrix
if degreeDistType == "uniform":
    k = simplexUtilities.generateUniformDegreeSequence(n, minDegree, maxDegree)
elif degreeDistType == "power-law":
    k = simplexUtilities.generatePowerLawDegreeSequence(n, minDegree, maxDegree, r)
elif degreeDistType == "poisson":
    k = simplexUtilities.generatePoissonDegreeSequence(n, meanDegree)

plt.figure()
plt.hist(k, bins=50)
plt.show()
A = simplexUtilities.generateConfigModelAdjacency(k)

# Calculate values needed in critical value calculation
meanDegree = simplexUtilities.meanPowerOfDegree(k, 1)
meanSquaredDegree =  simplexUtilities.meanPowerOfDegree(k, 2)
meanCubedDegree =  simplexUtilities.meanPowerOfDegree(k, 3)

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
    betaCrit = meanDegree/meanSquaredDegree*gamma
    alphaCrit = meanCubedDegree/(meanDegree**3 * meanSimplexDegree)*gamma

else:
    [simplexList, simplexIndices] = simplexUtilities.generateConfigModelSimplexList(k, simplexSize)
    # epidemic parameters
    betaCrit = meanDegree/meanSquaredDegree*gamma
    alphaCrit = (meanDegree**2)*meanCubedDegree/(meanSquaredDegree**3)*gamma

print("beta critical is " + str(betaCrit))
print("alpha critical is " + str(alphaCrit))

beta = betaCritFraction*betaCrit

start = time.time()
averageInfection, endState = simplexContagion.microscopicSimplexSISDynamics(A, simplexList, simplexIndices, gamma, beta, alpha, x0, timesteps, dt)
end = time.time()
print('The elapsed time is ' + str(end-start) + 's')

plt.figure()
plt.plot(np.linspace(0, (timesteps-1)*dt, timesteps), averageInfection)
plt.xlabel("time")
plt.ylabel(r"$\langle I\rangle$")
plt.show()
