import simplexTheory
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
minDegree = 50
maxDegreeList = np.linspace(100, 1000, 37)
rList = np.linspace(2.5,4.0,31)

numProcesses = len(os.sched_getaffinity(0))
print("Number of cores is " + str(numProcesses))
digits = 4
tolerance = 0.0001
minAlpha = 0.0
maxAlpha = 0.1
m = np.size(maxDegreeList,0)
n = np.size(rList,0)

firstOrderAlphaCritGrid = np.zeros([m,n])

argList = list()
for i in range(m):
    for j in range(n):
        degreeHist = generateTheoreticalDegreeHist(minDegree, int(maxDegreeList[i]), type, r=rList[j])

        meanDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 1)
        meanSquaredDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 2)
        meanCubedDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 3)
        meanSimplexDegree = meanDegree
        if isDegreeCorrelated:
            firstOrderAlphaCritGrid[i,j] = meanCubedDegree*meanDegree**2/meanSquaredDegree**3*gamma
        else:
            firstOrderAlphaCritGrid[i,j] = meanCubedDegree/meanDegree**4*gamma


        minBeta = 0.5*meanDegree/meanSquaredDegree*gamma
        maxBeta = 1.5*meanDegree/meanSquaredDegree*gamma

        argList.append((gamma, minBeta, maxBeta, minAlpha, maxAlpha, degreeHist, meanSimplexDegree, isDegreeCorrelated, digits, tolerance))

start = time.time()
with mp.Pool(processes=numProcesses) as pool:
    alphaCritList = pool.starmap(calculateTheoreticalCriticalAlpha, argList)
print("Time elapsed is " + str(time.time() - start) + "s", flush=True)

alphaCritGrid = np.reshape(alphaCritList, (m,n))

xMin = np.min(rList)
xMax = np.max(rList)
yMin = np.min(maxDegreeList)
yMax = np.max(maxDegreeList)

with open('alphaCrit' + datetime.now().strftime("%m%d%Y-%H%M%S"), 'wb') as file:
    pickle.dump([xMin, xMax, yMin, yMax, alphaCritGrid, firstOrderAlphaCritGrid], file)
