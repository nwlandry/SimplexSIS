import simplexUtilities
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from datetime import *

filename = 'Non-Random Degree/equilibriaData_power-law_r=3_dep_nonrandom_degree'
outputFilename = 'hysteresis_r=3_dep_nonrandom_degree'
simulationLabel = "Power law, r=3, Dependent (Simulation)"
theoreticalLabel = "Power law, r=3, Dependent (Theory)"
#filename = 'Archive-Data/equilibriaData11112019-002636'
with open(filename, 'rb') as file:
    data = pickle.load(file)

gamma = data[0]
beta = data[1]
alpha = data[2]
equilibria = data[3]
degreeSequence = data[4] # This is the full list of equilibria if it's an ensemble run
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
numBetaPoints = 60
numAlphaPoints = 30

betaTheory = np.linspace(min(beta),max(beta), numBetaPoints)
alphaTheory = np.linspace(min(alpha),max(alpha), numAlphaPoints)
hysteresisTheory = list()
for a in alphaTheory:
    hysteresisTheory.append(visualizeData.calculateTheoreticalHysteresis(gamma, betaTheory, a, minDegree, maxDegree, meanSimplexDegree, degreeSequence=degreeSequence, isIndependent=isIndependent, type="power-law", r=r, option="infinity", digits=4))

hysteresisSim = list()
for i in range(len(alpha)):
    hysteresisSim.append(visualizeData.calculateHysteresis(equilibria[i], beta, option='infinity'))

if not isinstance(degreeSequence[0], list):
    meanDegree = simplexUtilities.meanPowerOfDegree(degreeSequence, 1)
    meanSquaredDegree = simplexUtilities.meanPowerOfDegree(degreeSequence, 2)
    meanCubedDegree = simplexUtilities.meanPowerOfDegree(degreeSequence, 3)

if isIndependent:
    alphaCrit = meanSquaredDegree/(meanDegree**2*meanSimplexDegree)*gamma
else:
    alphaCrit = meanCubedDegree*meanDegree**2/(meanSquaredDegree**3)*gamma

with open(outputFilename, 'wb') as file:
    pickle.dump([alpha, hysteresisSim, alphaTheory, hysteresisTheory, simulationLabel, theoreticalLabel, alphaCrit], file)

plt.figure()
plt.plot(alpha, hysteresisTheory, 'k', label="Theory")
plt.plot(alpha, hysteresisSim, 'bo-', label="Simulation")
plt.scatter(alphaCrit, 0, s=50, color='red')
plt.xlabel(r"$\alpha$")
plt.ylabel("Hysteresis (sup-norm)")
plt.legend()
plt.plot()
plt.show()
