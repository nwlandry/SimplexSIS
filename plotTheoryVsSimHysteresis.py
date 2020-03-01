import simplexTheory
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *


#filename = 'Poster/uniform_indep'
#filename = 'equilibriaData12262019-000110'
#filename = 'Poster/power-law_r=4_indep'
filename = 'equilibriaData_power-law_r=3_indep_final'
filename = 'equilibriaData02292020-015156'
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
if len(degreeSequence[0]) != 1: # set degree sequence to none if "degree"
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

betaTheory = np.linspace(min(beta),max(beta), 41)

hysteresisTheory = list()
for i in range(len(alpha)):
    hysteresisTheory.append(visualizeData.calculateTheoreticalHysteresis(gamma, betaTheory, alpha[i], minDegree, maxDegree, meanSimplexDegree, degreeSequence=degreeSequence, isIndependent=isIndependent, type="power-law", r=r, option="infinity", digits=4))

hysteresisSim = list()
for i in range(len(alpha)):
    hysteresisSim.append(visualizeData.calculateHysteresis(equilibria[i], beta, option='infinity'))

plt.figure()
plt.plot(alpha, hysteresisTheory, 'k', label="Theory")
plt.plot(alpha, hysteresisSim, 'bo-', label="Simulation")
plt.xlabel(r"$\alpha$")
plt.ylabel("Hysteresis (sup-norm)")
plt.legend()
plt.plot()
plt.show()
