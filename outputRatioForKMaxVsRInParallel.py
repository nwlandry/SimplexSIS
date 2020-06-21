import simplexTheory
import simplexVisualize
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *
import simplexUtilities
import multiprocessing as mp
from datetime import datetime
import time
import os

gamma = 2
isDegreeCorrelated = False
type = "power-law"
minDegree = 50
maxDegreeList = np.linspace(100, 1000, 37)
exponentList = np.linspace(2.5,4.0,31)

numProcesses = len(os.sched_getaffinity(0))
print("Number of cores is " + str(numProcesses))
digits = 4
tolerance = 0.0001
option = "fast"
minAlpha = 0.0
maxAlpha = 0.1
m = np.size(maxDegreeList,0)
n = np.size(exponentList,0)

betaCritGrid = np.zeros([m,n])
expansionRatioGrid = np.zeros([m,n])

argList = list()
for i in range(m):
    for j in range(n):
        degreeHist = generateTheoreticalDegreeHist(minDegree, int(maxDegreeList[i]), type, exponent=exponentList[j])

        meanDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 1)
        meanSquaredDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 2)
        meanSimplexDegree = meanDegree

        betaCritGrid[i,j] = meanDegree/meanSquaredDegree*gamma

        expansionRatioGrid[i,j] = simplexUtilities.criticalBeta3(degreeHist, isDegreeCorrelated, gamma)/betaCritGrid[i,j]

        argList.append((gamma, betaCritGrid[i,j], minAlpha, maxAlpha, degreeHist, meanSimplexDegree, isDegreeCorrelated, digits, tolerance, option))

start = time.time()
with mp.Pool(processes=numProcesses) as pool:
    alphaCritList = pool.starmap(calculateTheoreticalCriticalAlpha, argList)
print("Time elapsed is " + str(time.time() - start) + "s", flush=True)

alphaCritGrid = np.reshape(alphaCritList, (m,n))
ratioGrid = np.divide(alphaCritGrid,betaCritGrid)

xMin = np.min(exponentList)
xMax = np.max(exponentList)
yMin = np.min(maxDegreeList)
yMax = np.max(maxDegreeList)

with open('ratio' + datetime.now().strftime("%m%d%Y-%H%M%S"), 'wb') as file:
    pickle.dump([xMin, xMax, yMin, yMax, ratioGrid, expansionRatioGrid], file)
