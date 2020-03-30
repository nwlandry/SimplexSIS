import simplexTheory
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *
import multiprocessing as mp
from datetime import datetime
import time

gamma = 2
isIndependent = False
type = "power-law"
minDegree = 50
maxDegree = 1000
gamma = 3

numSimulations = 1000
numProcesses = 24
digits = 4
tolerance = 0.0001
numBetaPoints = 50
minAlpha = 0.0
maxAlpha = 0.1

argList = list()
for i in range(numSimulations):
    degreeSequence = generatePowerLawDegreeSequence(numPoints, minDegree, maxDegree, r, isRandom=True)
    degreeHist = degreeSequenceToHist(degreeSequence)

    meanDegree = sum([k*prob for k, prob in degreeHist])
    meanSquaredDegree = sum([k**2*prob for k, prob in degreeHist])
    meanSimplexDegree = meanDegree

    betaTheory = np.linspace(0.5*meanDegree/meanSquaredDegree*gamma, 1.5*meanDegree/meanSquaredDegree*gamma, numBetaPoints)

    argList.append((gamma, betaTheory, minAlpha, maxAlpha, degreeHist, meanSimplexDegree, isIndependent, "infinity", digits, tolerance))

with mp.Pool(processes=numProcesses) as pool:
    alphaCritList = pool.starmap(calculateTheoreticalCriticalAlpha, argList)

with open('alphaCritList' + datetime.now().strftime("%m%d%Y-%H%M%S"), 'wb') as file:
    pickle.dump(alphaCritList, file)
