import simplexUtilities
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from datetime import *
from simplexTheory import *

gamma = 2
minDegree = 50
maxDegree = 450
isIndependent = True
type = "power-law"

r = 3.0
digits = 5
tolerance = 0.0001
minAlpha = 0
maxAlpha = 0.6

degreeHist = generateTheoreticalDegreeHist(minDegree, maxDegree, type, r=r)

meanDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 1)
meanSquaredDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 2)
meanSimplexDegree = meanDegree

minBeta = 0.5*meanDegree/meanSquaredDegree*gamma
maxBeta = 1.5*meanDegree/meanSquaredDegree*gamma
alphaCrit = calculateTheoreticalCriticalAlpha(gamma, minBeta, maxBeta, minAlpha, maxAlpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, option="infinity", digits=digits, tolerance=tolerance)
print(alphaCrit)
