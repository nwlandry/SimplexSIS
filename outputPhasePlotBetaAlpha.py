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

minBeta = 0.5*meanDegree/meanSquaredDegree*gamma
maxBeta = 1.5*meanDegree/meanSquaredDegree*gamma
beta = np.linspace(minBeta, maxBeta, 51)

minAlpha = 0.0
maxAlpha = 0.1
alpha = np.linspace(minAlpha, maxAlpha, 51)

numProcesses = mp.cpu_count()
digits = 4
tolerance = 0.0001

m = np.size(alpha, 0)
n = np.size(beta, 0)

argList = list()

for i in range(m):
    for j in range(n):
        argList.append((gamma, beta[j], alpha[i], degreeHist, meanSimplexDegree, isIndependent, majorityVote, digits))

with mp.Pool(processes=numProcesses) as pool:
    phaseGridList = pool.starmap(getPhase, argList)

phaseGrid = np.reshape(phaseGridList, (m,n))

xMin = np.min(beta)
xMax = np.max(beta)
yMin = np.min(alpha)
yMax = np.max(alpha)

with open('phasePlot' + datetime.now().strftime("%m%d%Y-%H%M%S"), 'wb') as file:
    pickle.dump([xMin, xMax, yMin, yMax, phaseGrid], file)
