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

meanDegree = sum([k*prob for k, prob in degreeHist])
meanSquaredDegree = sum([k**2*prob for k, prob in degreeHist])
meanCubedDegree = sum([k**3*prob for k, prob in degreeHist])
meanSimplexDegree = meanDegree

betaTheory = np.linspace(0.5*meanDegree/meanSquaredDegree*gamma, 1.5*meanDegree/meanSquaredDegree*gamma, 9)
alphaCrit = calculateTheoreticalCriticalAlpha(gamma, betaTheory, minAlpha, maxAlpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, option="infinity", digits=digits, tolerance=tolerance)
print(alphaCrit)
