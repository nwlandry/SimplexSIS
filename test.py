import simplexTheory
from simplexTheory import getPhase
import simplexVisualize
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *
import multiprocessing as mp
from datetime import datetime
import time

gamma = 2
isIndependent = True
type = "power-law"
minDegree = 50
maxDegree = 100
r = 4.0
majorityVote = True

degreeHist = generateTheoreticalDegreeHist(minDegree, maxDegree, type, r=r)

meanDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 1)
meanSquaredDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 2)
meanSimplexDegree = meanDegree

minBeta = 0*meanDegree/meanSquaredDegree*gamma
maxBeta = 1.5*meanDegree/meanSquaredDegree*gamma
beta = np.linspace(minBeta, maxBeta, 10)

minAlpha = 0
maxAlpha = 0.1
alpha = np.linspace(minAlpha, maxAlpha, 10)

numProcesses = mp.cpu_count()
print(numProcesses)
digits = 4
tolerance = 0.0001

m = np.size(alpha, 0)
n = np.size(beta, 0)

phaseGrid = np.zeros([m,n])

healing = True
for i in range(m):
    for j in range(n):
        phaseGrid[i,j] = simplexTheory.getPhase(gamma, beta[j], alpha[i], degreeHist, meanSimplexDegree, isIndependent, majorityVote, healing, digits)

xMin = np.min(beta)
xMax = np.max(beta)
yMin = np.min(alpha)
yMax = np.max(alpha)

plt.figure()
c = plt.imshow(np.flipud(phaseGrid), interpolation="none", cmap="Reds", extent=[xMin, xMax, yMin, yMax], aspect="auto", vmin=1, vmax=3)
plt.xlabel(r"$\beta_2$")
plt.ylabel(r"$\beta_3$")
plt.plot()
plt.show()
