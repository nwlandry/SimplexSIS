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
# graph parameters
k0 = 50
r = 4 # power law exponent
minDeg = 5
maxDeg = 10
n = 11
#nSimplices = 1000
simplexSize = 3

# Epidemic parameters
initialFraction = 0.5
x0 = np.random.choice([0, 1], size=n, p=[1-initialFraction, initialFraction])

#simulation parameters
timesteps = 1000
dt = 0.1
avgLength = int(0.4*timesteps)
numBetaPts = 4
numAlphaPts = 2

# generate adjacency matrices and simplex list
k = simplexUtilities.generatePowerLawDegreeSequence(numPoints, k0, n, r)
A = simplexUtilities.generateConfigModelAdjacency(k)
[simplexList, simplexIndices] = simplexUtilities.generateConfigModelSimplexList(k, simplexSize)

# Calculate critical values
kAvg = A.sum(axis=0).mean(axis=1).item(0,0)
kSquaredAvg = ((A.sum(axis=0)*A.sum(axis=1)).mean(axis=0).mean(axis=1)/n).item(0,0)
kCubedAvg = ((A.sum(axis=1)*A.sum(axis=0)*A.sum(axis=1)).mean(axis=0)/n).item(0,0)

print(kAvg)
print(kSquaredAvg)
print(kCubedAvg)

# epidemic parameters
gamma = 2
betaCrit = kAvg/kSquaredAvg*gamma
alphaCrit = (kAvg**2)*kCubedAvg/(kSquaredAvg**3)*gamma

print(betaCrit)
print(alphaCrit)

beta = np.concatenate([np.linspace(0, 2*betaCrit, numBetaPts),
                      np.linspace(2*betaCrit*(1-1.0/(numBetaPts-1)), 0, numBetaPts-1)])
alpha = np.linspace(0, 2*alphaCrit, numAlphaPts)

start = time.time()
equilibria = simplexContagion.generateSISEquilibria(A, simplexList, simplexIndices, gamma, beta, alpha, x0, timesteps, dt, avgLength)
end = time.time()
print('The elapsed time is ' + str(end-start) + 's')

with open('equilibriaData' + datetime.now().strftime("%m%d%Y-%H%M%S"), 'wb') as file:
    pickle.dump(equilibria, file)
    pickle.dump(alpha, file)
    pickle.dump(beta, file)
