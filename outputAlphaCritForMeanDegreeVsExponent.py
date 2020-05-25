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
maxDegree = 500
meanDegreeList = np.linspace(50, 200, 16)
exponentList = np.linspace(2.501,4.001,31)

numProcesses = len(os.sched_getaffinity(0))
print("Number of cores is " + str(numProcesses))
digits = 4
tolerance = 0.0001
option = "fast"
minAlpha = 0.0
maxAlpha = 0.4
initialGuess = 10
m = np.size(meanDegreeList,0)
n = np.size(exponentList,0)

firstOrderAlphaCritGrid = np.zeros([m,n])

argList = list()
for i in range(m):
    for j in range(n):
        minDegree = fsolve(simplexTheory.avgOfPowerLawEqn, initialGuess, args=(maxDegree, exponentList[j], meanDegreeList[i]))
        degreeHist = generateTheoreticalDegreeHist(int(minDegree), maxDegree, type, power=exponentList[j])

        meanDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 1)
        meanSquaredDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 2)
        meanCubedDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 3)
        meanSimplexDegree = meanDegree

        if isDegreeCorrelated:
            firstOrderAlphaCritGrid[i,j] = meanCubedDegree*meanDegree**2/meanSquaredDegree**3*gamma
        else:
            firstOrderAlphaCritGrid[i,j] = meanCubedDegree/meanDegree**4*gamma

        betaCrit = meanDegree/meanSquaredDegree*gamma

        argList.append((gamma, betaCrit, minAlpha, maxAlpha, degreeHist, meanSimplexDegree, isDegreeCorrelated, digits, tolerance, option))

start = time.time()
with mp.Pool(processes=numProcesses) as pool:
    alphaCritList = pool.starmap(calculateTheoreticalCriticalAlpha, argList)
print("Time elapsed is " + str(time.time() - start) + "s", flush=True)

alphaCritGrid = np.reshape(alphaCritList, (m,n))

xMin = np.min(exponentList)
xMax = np.max(exponentList)
yMin = np.min(meanDegreeList)
yMax = np.max(meanDegreeList)

with open('alphaCrit' + datetime.now().strftime("%m%d%Y-%H%M%S"), 'wb') as file:
    pickle.dump([xMin, xMax, yMin, yMax, alphaCritGrid, firstOrderAlphaCritGrid], file)
