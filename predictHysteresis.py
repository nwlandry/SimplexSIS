import simplexTheory
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *

gamma = 2
minBeta = 0
maxBeta = 0.05
numBetaPoints = 30
minAlpha = 0.04
maxAlpha = 0.05
numAlphaPoints = 20
beta = np.linspace(minBeta, maxBeta, numBetaPoints)
alpha = np.linspace(minAlpha, maxAlpha, numAlphaPoints)

minDegree = 50
maxDegree = 450
isIndependent = True
type = "power-law"
meanSimplexDegree = 100
r = 3.0
digits = 4

hysteresisTheory = list()
for i in range(len(alpha)):
    hysteresisTheory.append(visualizeData.calculateTheoreticalHysteresis(gamma, beta, alpha[i], minDegree, maxDegree, meanSimplexDegree, degreeSequence=None, isIndependent=isIndependent, type=type, r=r, option="infinity", digits=digits))

plt.figure()
plt.plot(alpha, hysteresisTheory, 'k')
plt.xlabel(r"$\alpha$")
plt.ylabel("Hysteresis (sup-norm)")
plt.plot()
plt.show()
