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
filename = 'equilibriaData_power-law_r=4_dep_final'
#filename = 'Archive-Data/equilibriaData10252019-225011'
#filename = 'Archive-Data/equilibriaData11112019-002636'
with open(filename, 'rb') as file:
    data = pickle.load(file)

gamma = data[0]
beta = data[1]
alpha = data[2]
equilibria = data[3]
betaCrit = data[4]
alphaCrit = data[5]
meanDegree = data[6]
meanSquaredDegree = data[7]
meanCubedDegree = data[8]
meanSimplexDegree = data[9]
degreeSequence = data[10]

#degreeSequence = None
minDegree = 67
maxDegree = 450
isIndependent = False
type = "power-law"
meanSimplexDegree = 100
r = 4.0
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
