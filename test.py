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

degreeHist = generateTheoreticalDegreeHist(minDegree, maxDegree, type, r=r)

meanDegree = sum([k*prob for k, prob in degreeHist])
meanSquaredDegree = sum([k**2*prob for k, prob in degreeHist])
meanCubedDegree = sum([k**3*prob for k, prob in degreeHist])
meanSimplexDegree = meanDegree
print(meanDegree)
print(meanSquaredDegree)
print(meanCubedDegree)
print(sum([prob for k,prob in degreeHist]))

minAlpha = 0
maxAlpha = 0.6
betaTheory = np.linspace(0.5*meanDegree/meanSquaredDegree*gamma, 1.5*meanDegree/meanSquaredDegree*gamma, 9)
alphaCrit = calculateTheoreticalCriticalAlpha(gamma, betaTheory, minAlpha, maxAlpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, option="infinity", digits=4, tolerance=0.001)
print(alphaCrit)
if isIndependent:
    alphaCrit2 = 4*(meanSquaredDegree*gamma-meanDegree**2)/(meanSquaredDegree*meanSimplexDegree)
else:
    alphaCrit2 = meanCubedDegree*meanDegree**2/(meanSquaredDegree**3)*gamma

print(alphaCrit2)

# fractionAlphaCrit = 1
#
# alpha = fractionAlphaCrit*alphaCrit
# print(alphaCrit)
# print(Stop)
# beta = np.linspace(0.5*meanDegree/meanSquaredDegree*gamma,1.5*meanDegree/meanSquaredDegree*gamma,200)
#
# visualizeData.plotTheoreticalInfectionCurves(gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, degreeHist=degreeHist, isIndependent=isIndependent, type=type, r=r, digits=4)
