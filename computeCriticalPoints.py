import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs
import simplexUtilities
import simplexContagion

# graph parameters
k0 = 50
r = 2 # power law exponent
minDeg = 50
maxDeg = 100
n = 1000
simplexSize = 3
isIndependentUniform = False
degreeDistType = "uniform"
gamma = 2

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

alphaCritIndep = kCubedAvg/kAvg**4*gamma
alphaCritDependent = kAvg**2*kCubedAvg/kSquaredAvg**3*gamma
print("The alpha critical for independent uniformly chosen simplices is " + str(alphaCritIndep))
print("The alpha critical for simplices correlated to degree is " + str(alphaCritDependent))
