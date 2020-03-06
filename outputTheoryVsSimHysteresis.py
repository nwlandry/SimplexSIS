import simplexTheory
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *

filename = 'equilibriaData_power-law_r=4_dep_final'
outputFilename = 'hysteresis_r=4_dep_final'
simulationLabel = "Power law, r=4, Dependent (Simulation)"
theoreticalLabel = "Power law, r=4, Dependent (Theory)"
#filename = 'Archive-Data/equilibriaData11112019-002636'
with open(filename, 'rb') as file:
    data = pickle.load(file)

gamma = data[0]
beta = data[1]
alpha = data[2]
equilibria = data[3]
#degreeSequence = data[4] # This is the full list of equilibria if it's an ensemble run
degreeSequence = data[10]
#isIndependent = data[5]
#type = data[6]
#r = data[7]
isIndependent = False
r = 4
type = "power-law"
if isinstance(degreeSequence[0], list) : # set degree sequence to none if "degree"
    degreeSequence = None
    meanDegree = data[8]
    meanSquaredDegree = data[9]
    meanCubedDegree = data[10]
    meanSimplexDegree = data[11]
    minDegree = data[12]
    maxDegree = data[13]
else:
    meanSimplexDegree = data[9]
    #meanSimplexDegree = data[8]
    minDegree = min(degreeSequence)
    maxDegree = max(degreeSequence)

digits = 5

betaTheory = np.linspace(min(beta),max(beta), 61)

hysteresisTheory = list()
for i in range(len(alpha)):
    hysteresisTheory.append(visualizeData.calculateTheoreticalHysteresis(gamma, betaTheory, alpha[i], minDegree, maxDegree, meanSimplexDegree, degreeSequence=degreeSequence, isIndependent=isIndependent, type="power-law", r=r, option="infinity", digits=4))

hysteresisSim = list()
for i in range(len(alpha)):
    hysteresisSim.append(visualizeData.calculateHysteresis(equilibria[i], beta, option='infinity'))

with open(outputFilename, 'wb') as file:
    pickle.dump([alpha, hysteresisSim, hysteresisTheory, simulationLabel, theoreticalLabel], file)

plt.figure()
plt.plot(alpha, hysteresisTheory, 'k', label="Theory")
plt.plot(alpha, hysteresisSim, 'bo-', label="Simulation")
plt.xlabel(r"$\alpha$")
plt.ylabel("Hysteresis (sup-norm)")
plt.legend()
plt.plot()
plt.show()
