import simplexUtilities
import simplexTheory
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from datetime import *

filename = 'Non-Random Degree/equilibriaData_uniform_indep_nonrandom_degree'
outputFilename = 'hysteresis_uniform_indep_nonrandom_degree'
simulationLabel = "Uniform, Independent (Simulation)"
theoreticalLabel = "Uniform, Independent (Theory)"
#filename = 'Archive-Data/equilibriaData11112019-002636'
with open(filename, 'rb') as file:
    data = pickle.load(file)

gamma = data[0]
beta = data[1]
alphaSim = data[2]
equilibria = data[3]
degreeSequence = data[4] # This is the full list of equilibria if it's an ensemble run
# degreeSequence = data[10]

isIndependent = data[5]
type = data[6]
r = data[7]
if isinstance(degreeSequence[0], list) : # set degree sequence to none if "degree"
    degreeSequence = None
    meanDegree = data[8]
    meanSquaredDegree = data[9]
    meanCubedDegree = data[10]
    meanSimplexDegree = data[11]
    minDegree = data[12]
    maxDegree = data[13]
else:
    meanSimplexDegree = data[8]
    minDegree = min(degreeSequence)
    maxDegree = max(degreeSequence)

digits = 5
tolerance = 0.0001
numAlphaPoints = 50

meanDegree = simplexUtilities.meanPowerOfDegree(degreeSequence, 1)
meanSquaredDegree = simplexUtilities.meanPowerOfDegree(degreeSequence, 2)
#meanSimplexDegree = meanDegree
degreeHist = simplexTheory.degreeSequenceToHist(degreeSequence)

betaCrit = meanDegree/meanSquaredDegree*gamma
minBeta = 0.5*betaCrit
maxBeta = 1.5*betaCrit

alphaTheory = np.linspace(min(alpha),max(alpha), numAlphaPoints)
hysteresisTheory = list()
for a in alphaTheory:
    hysteresisTheory.append(simplexTheory.calculateTheoreticalHysteresis(gamma, minBeta, maxBeta, a, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, option="infinity", digits=digits, tolerance=tolerance))

hysteresisSim = list()
for i in range(len(alphaSim)):
    hysteresisSim.append(visualizeData.calculateHysteresis(equilibria[i], beta, option='infinity'))

if not isinstance(degreeSequence[0], list):
    meanDegree = simplexUtilities.meanPowerOfDegree(degreeSequence, 1)
    meanSquaredDegree = simplexUtilities.meanPowerOfDegree(degreeSequence, 2)
    meanCubedDegree = simplexUtilities.meanPowerOfDegree(degreeSequence, 3)

    print(hysteresisTheory)

with open(outputFilename, 'wb') as file:
    pickle.dump([alphaSim, hysteresisSim, alphaTheory, hysteresisTheory, simulationLabel, theoreticalLabel], file)

plt.figure()
plt.plot(alphaTheory, hysteresisTheory, 'k', label="Theory")
plt.plot(alphaSim, hysteresisSim, 'bo-', label="Simulation")
plt.xlabel(r"$\alpha$")
plt.ylabel("Hysteresis (sup-norm)")
plt.legend()
plt.plot()
plt.show()
