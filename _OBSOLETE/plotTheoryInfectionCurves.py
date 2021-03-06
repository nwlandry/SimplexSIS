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
isDegreeCorrelated = True
type = "power-law"
minDegree = 50
maxDegree = 450
exponent = 4.0
digits = 4
tolerance = 0.0001
degreeHist = generateTheoreticalDegreeHist(minDegree, maxDegree, type, exponent=exponent)
meanDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 1)
meanSquaredDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 2)
meanSimplexDegree = meanDegree

minAlpha = 0
maxAlpha = 0.06
betaCrit = meanDegree/meanSquaredDegree*gamma
minBeta = 0.5*betaCrit
maxBeta = 1.5*betaCrit

numBetaPoints = 50

#alphaCrit = calculateTheoreticalCriticalAlpha(gamma, betaCrit, minAlpha, maxAlpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits, tolerance=tolerance)
alphaCrit = 0.1
alphaCritFraction = [1.0]

beta = np.linspace(minBeta, maxBeta, numBetaPoints)
simplexVisualize.plotTheoreticalInfectionCurves(gamma, beta, alphaCritFraction, alphaCrit, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits)
