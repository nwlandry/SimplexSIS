import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs
import random
import simplexUtilities
import simplexContagion
import simplexTheory
import pickle
from datetime import datetime
import time
import multiprocessing as mp

# graph parameters
k0 = 50
r = 4 # power law exponent
minDeg = 50
maxDeg = 100
n = 1000
simplexSize = 3
isIndependentUniform = False
degreeDistType = "power-law"

# Epidemic parameters
initialFraction = 0.5
x0 = np.random.choice([0, 1], size=n, p=[1-initialFraction, initialFraction])

#simulation parameters
timesteps = 10000
dt = 0.1
avgLength = int(0.4*timesteps)
numBetaPts = 31
numAlphaPts = 8
startAlphaCritFraction = 0
endAlphaCritFraction = 4
startBetaCritFraction = 0.5
endBetaCritFraction = 1.5
numProcesses = numAlphaPts


# generate degree sequence and adjacency matrix
if degreeDistType == "uniform":
    k = simplexUtilities.generateUniformDegreeSequence(n, minDeg, maxDeg)
elif degreeDistType == "power-law":
    k = simplexUtilities.generatePowerLawDegreeSequence(n, k0, n, r)

degreeHist = simplexTheory.degreeSequenceToHist(k)
print(degreeHist)
sum = 0
for degree in degreeHist:
    sum = sum +degree[0]*degree[1]
print(len(degreeHist))
print(sum)
