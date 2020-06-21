import simplexUtilities
import simplexTheory
import simplexVisualize
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from datetime import *

filename = 'Paper Figures/equilibriaData_r=3_degreecorrelated'
filename = 'equilibriaData06142020-212240'
outputFilename = 'bistability_r=3_degreecorrelated'
simulationLabel = r"$P(k)\propto k^{-3}$ (Simulation)"
theoreticalLabel = r"$P(k)\propto k^{-3}$ (Theory)"
#filename = 'Archive-Data/equilibriaData11112019-002636'
with open(filename, 'rb') as file:
    data = pickle.load(file)

gamma = data[0]
beta = data[1]
alphaSim = data[2]
equilibria = data[3]
degreeSequence = data[4]
isDegreeCorrelated = data[5]
type = data[6]
exponent = data[7]
meanSimplexDegree = data[8]
minDegree = min(degreeSequence)
maxDegree = max(degreeSequence)

digits = 5
tolerance = 0.0001
numAlphaPoints = 100

meanDegree = simplexUtilities.meanPowerOfDegree(degreeSequence, 1)
meanSquaredDegree = simplexUtilities.meanPowerOfDegree(degreeSequence, 2)
degreeHist = simplexTheory.degreeSequenceToHist(degreeSequence)

betaCrit = meanDegree/meanSquaredDegree*gamma

alphaTheory = np.linspace(min(alphaSim),max(alphaSim), numAlphaPoints)
bistabilityTheory = list()
for a in alphaTheory:
    bistabilityTheory.append(simplexTheory.calculateTheoreticalBistability(gamma, betaCrit, a, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits, tolerance=tolerance))

bistabilitySim = list()
for i in range(len(alphaSim)):
    bistabilitySim.append(simplexVisualize.calculateBistability(equilibria[i], beta, option='infinity'))

with open(outputFilename, 'wb') as file:
    pickle.dump([alphaSim, bistabilitySim, alphaTheory, bistabilityTheory, simulationLabel, theoreticalLabel], file)

plt.figure()
plt.plot(alphaTheory, bistabilityTheory, 'k', label="Theory")
plt.plot(alphaSim, bistabilitySim, 'bo-', label="Simulation")
plt.xlabel(r"$\beta_3$")
plt.ylabel(r"$B(\beta_3)$")
plt.legend()
plt.plot()
plt.show()
