import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs
import simplexUtilities
import simplexContagion
import matplotlib.pyplot as plt
# graph parameters
k0 = 10
r = 4 # power law exponent
minDeg = 50
maxDeg = 150
n = 1000
degreeDistType = "uniform"
gamma = 2
numSims = 10000
alphaCritIndependent = list()
alphaCritDependent = list()
for i in range(numSims):
    # generate degree sequence and adjacency matrix
    if degreeDistType == "uniform":
        k = simplexUtilities.generateUniformDegreeSequence(n, minDeg, maxDeg)
    elif degreeDistType == "power-law":
        k = simplexUtilities.generatePowerLawDegreeSequence(n, k0, n, r)

    # Calculate critical values
    meanDegree = simplexUtilities.meanPowerOfDegree(k, 1)
    meanSquaredDegree =  simplexUtilities.meanPowerOfDegree(k, 2)
    meanCubedDegree =  simplexUtilities.meanPowerOfDegree(k, 3)
    alphaCritIndependent.append(meanCubedDegree/(meanDegree**4)*gamma)
    alphaCritDependent.append((meanDegree**2)*meanCubedDegree/(meanSquaredDegree**3)*gamma)

plt.figure()
plt.hist(alphaCritIndependent, bins=100)
plt.xlabel(r'$\alpha_{crit}$')
plt.ylabel("Instances")
plt.show()

plt.figure()
plt.hist(alphaCritDependent, bins=100)
plt.xlabel(r'$\alpha_{crit}$')
plt.ylabel("Instances")
plt.show()
