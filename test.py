import simplexUtilities
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from datetime import *
from simplexTheory import *
import time


gamma = 2
minDegree = 50
maxDegree = 450
isIndependent = True
type = "power-law"

r = 3.0
digits = 5

degreeHist = generateTheoreticalDegreeHist(minDegree, maxDegree, type, r=r)

meanDegree = sum([k*prob for k, prob in degreeHist])
meanSquaredDegree = sum([k**2*prob for k, prob in degreeHist])
meanCubedDegree = sum([k**3*prob for k, prob in degreeHist])
meanSimplexDegree = meanDegree
print(meanDegree)
print(meanSquaredDegree)
print(meanCubedDegree)

minAlpha = 0
maxAlpha = 0.06
minBeta = 0
maxBeta = 2*meanDegree/meanSquaredDegree*gamma
betaTheory = np.linspace(0.5*meanDegree/meanSquaredDegree*gamma, 1.5*meanDegree/meanSquaredDegree*gamma, 30)
#alphaCrit = calculateTheoreticalCriticalAlpha(gamma, betaTheory, minAlpha, maxAlpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, option="infinity", digits=4, tolerance=0.001)
#print(alphaCrit)
#alphaCrit2 = calculateTheoreticalCriticalAlphaFast(gamma, minBeta, maxBeta, minAlpha, maxAlpha, degreeHist, meanSimplexDegree=None, isIndependent=False, option="infinity", digits=4, tolerance=0.0001)
alpha = 0.0362

start = time.time()
orig = calculateTheoreticalHysteresisOriginal(gamma, betaTheory, alpha, degreeHist, meanSimplexDegree=None, isIndependent=False, option="infinity", digits=4)
end = time.time()
print('The elapsed time is ' + str(end-start) + 's')
print(orig)

start = time.time()
test = calculateTheoreticalHysteresis(gamma, minBeta, maxBeta, alpha, degreeHist, meanSimplexDegree=None, isIndependent=False, option="infinity", digits=4, tolerance=0.0001)
end = time.time()
print('The elapsed time is ' + str(end-start) + 's')
print(test)

# fractionAlphaCrit = 1
#
# alpha = fractionAlphaCrit*alphaCrit
# print(alphaCrit)
# print(Stop)
# beta = np.linspace(0.5*meanDegree/meanSquaredDegree*gamma,1.5*meanDegree/meanSquaredDegree*gamma,200)
#
# visualizeData.plotTheoreticalInfectionCurves(gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, degreeHist=degreeHist, isIndependent=isIndependent, type=type, r=r, digits=4)
