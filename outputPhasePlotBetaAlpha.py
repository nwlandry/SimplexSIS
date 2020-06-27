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
import os

gamma = 2
isDegreeCorrelated = False
type = "power-law"
minDegree = 67
maxDegree = 1000
exponent = 4.0
majorityVote = True
healing = False

degreeHist = generateTheoreticalDegreeHist(minDegree, maxDegree, type, exponent=exponent)

meanDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 1)
meanSquaredDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 2)
meanSimplexDegree = meanDegree

minBeta = 0.5*meanDegree/meanSquaredDegree*gamma
maxBeta = 1.5*meanDegree/meanSquaredDegree*gamma
beta = np.linspace(minBeta, maxBeta, 101)

minAlpha = 0.0
maxAlpha = 0.06
alpha = np.linspace(minAlpha, maxAlpha, 101)

numProcesses = len(os.sched_getaffinity(0))
digits = 4
tolerance = 0.0001

m = np.size(alpha, 0)
n = np.size(beta, 0)

argList = list()

for i in range(m):
    for j in range(n):
        argList.append((gamma, beta[j], alpha[i], degreeHist, meanSimplexDegree, isDegreeCorrelated, majorityVote, healing, digits))

with mp.Pool(processes=numProcesses) as pool:
    phaseGridList = pool.starmap(getPhase, argList)

phaseGrid = np.reshape(phaseGridList, (m,n))

xMin = np.min(beta)
xMax = np.max(beta)
yMin = np.min(alpha)
yMax = np.max(alpha)

with open('phasePlot' + datetime.now().strftime("%m%d%Y-%H%M%S"), 'wb') as file:
    pickle.dump([xMin, xMax, yMin, yMax, phaseGrid], file)
