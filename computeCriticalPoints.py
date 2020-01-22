import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs
import simplexUtilities
import simplexContagion

# graph parameters
k0 = 10
r = 4 # power law exponent
minDeg = 10
maxDeg = 1000
n = 10000
degreeDistType = "power-law"
gamma = 2

# generate degree sequence and adjacency matrix
if degreeDistType == "uniform":
    k = simplexUtilities.generateUniformDegreeSequence(n, minDeg, maxDeg)
elif degreeDistType == "power-law":
    k = simplexUtilities.generatePowerLawDegreeSequence(n, k0, n, r)

# Calculate critical values
meanDegree = simplexUtilities.meanPowerOfDegree(k, 1)
meanSquaredDegree =  simplexUtilities.meanPowerOfDegree(k, 2)
meanCubedDegree =  simplexUtilities.meanPowerOfDegree(k, 3)
print(meanDegree)
print(meanSquaredDegree)
print(meanCubedDegree)
alphaCritIndep = meanCubedDegree/(meanDegree**4)*gamma
alphaCritDependent = (meanDegree**2)*meanCubedDegree/(meanSquaredDegree**3)*gamma
print("The alpha critical for independent uniformly chosen simplices is " + str(alphaCritIndep))
print("The alpha critical for simplices correlated to degree is " + str(alphaCritDependent))
