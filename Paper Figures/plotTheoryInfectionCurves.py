import sys
sys.path.insert(1, 'C:/Users/nicho/Documents/GitHub/SimplexSIS')

import simplexTheory
import simplexVisualize
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *
import math

gamma = 2
isIndependent = True
type = "power-law"
minDegree = 50
maxDegree = 450
r = 4.0
digits = 6
tolerance = 0.0001
degreeHist = generateTheoreticalDegreeHist(minDegree, maxDegree, type, r=r)
meanDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 1)
meanSquaredDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 2)
meanSimplexDegree = meanDegree

minAlpha = 0
maxAlpha = 0.06
minBeta = 0*meanDegree/meanSquaredDegree*gamma
maxBeta = 2*meanDegree/meanSquaredDegree*gamma
numBetaPoints = 50

#alphaCrit = calculateTheoreticalCriticalAlpha(gamma, minBeta, maxBeta, minAlpha, maxAlpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, option="infinity", digits=digits, tolerance=tolerance)
alphaCrit = 0.1
print(alphaCrit)
alphaCritFraction = [1.0]

beta = np.linspace(minBeta, maxBeta, numBetaPoints)
simplexVisualize.plotTheoreticalInfectionCurves(gamma, beta, alphaCritFraction, alphaCrit, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, digits=digits)
