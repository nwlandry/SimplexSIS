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
isIndependent = False
type = "power-law"
minDegree = 50
maxDegree = 60
r = 4.0
digits = 4
tolerance = 0.000001
degreeHist = generateTheoreticalDegreeHist(minDegree, maxDegree, type, r=r)
meanDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 1)
meanSquaredDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 2)
meanCubedDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 3)
meanSimplexDegree = meanDegree

alphaCrit = meanCubedDegree*meanDegree**2/meanSquaredDegree**3*gamma
print(alphaCrit)
alpha = 0.04
minBeta = 0.85*meanDegree/meanSquaredDegree*gamma
maxBeta = 1.15*meanDegree/meanSquaredDegree*gamma
numBetaPoints = 50

hys = calculateTheoreticalHysteresisVisually(gamma, minBeta, maxBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, option="infinity", digits=digits, tolerance=tolerance, stopAtBistability=False)
print(hys)
#
# beta = np.linspace(minBeta, maxBeta, numBetaPoints)
# simplexVisualize.plotTheoreticalInfectionCurves(gamma, beta, [1.0], alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, digits=digits)
