import simplexTheory
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *

gamma = 2
minDegree = 50
maxDegrees = [100,150,200,250,300,350,400,450]
isIndependent = True
type = "power-law"
r = 3.0
digits = 4
tolerance = 0.001
numBetaPoints = 30
minAlpha = 0.04
maxAlpha = 0.08

critAlphaTheory = list()
critAlphaLinear = list()

for maxDegree in maxDegrees:
    degreeHist = generateTheoreticalDegreeHist(minDegree, maxDegree, type, r=r)
    meanDegree = sum([k*prob for k, prob in degreeHist])
    meanSquaredDegree = sum([k**2*prob for k, prob in degreeHist])
    meanCubedDegree = sum([k**3*prob for k, prob in degreeHist])
    meanSimplexDegree = meanDegree
    beta = np.linspace(0.5*meanDegree/meanSquaredDegree*gamma,1.5*meanDegree/meanSquaredDegree*gamma,numBetaPoints)
    if isIndependent:
        alphaCrit = meanCubedDegree/(meanDegree**3*meanSimplexDegree)*gamma
        alphaCrit = 4*(meanSquaredDegree*gamma-meanDegree**2)/(meanSquaredDegree*meanSimplexDegree)
        critAlphaLinear.append(alphaCrit)
    else:
        alphaCrit = meanCubedDegree*meanDegree**2/(meanSquaredDegree**3)*gamma
        critAlphaLinear.append(alphaCrit)

    alphaCrit = calculateTheoreticalCriticalAlpha(gamma, beta, minAlpha, maxAlpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, option="infinity", digits=digits, tolerance=tolerance)
    critAlphaTheory.append(alphaCrit)

    print(maxDegree)

plt.figure()
plt.plot(maxDegrees, critAlphaTheory, 'k', label='Mean Field')
plt.plot(maxDegrees, critAlphaLinear, 'b', label='Linearized')
plt.xlabel("Maximum degree")
plt.ylabel(r"Critical $\alpha$")
plt.legend(loc='lower right')
plt.plot()
plt.show()
