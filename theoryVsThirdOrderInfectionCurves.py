import simplexTheory
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *
import csv

#filename = 'Poster/uniform_indep'
#filename = 'equilibriaData12262019-000110'
#filename = 'Poster/power-law_r=4_indep'
#filename = 'equilibriaData_power-law_r=4_final'
#filename = 'Archive-Data/equilibriaData10252019-225011'
#filename = 'Archive-Data/equilibriaData11112019-002636'
mathematicaFile = 'indepData.csv'
# with open(filename, 'rb') as file:
#     data = pickle.load(file)
#
# gamma = data[0]
# beta = data[1]
# alpha = data[2]
# equilibria = data[3]



gamma = 2
minDegree = 67
maxDegree = 450
isIndependent = True
type = "power-law"
meanSimplexDegree = 100
r = 4.0
digits = 5

degreeHist = generateTheoreticalDegreeHist(minDegree, maxDegree, type, r=r)

index = 0

alpha = 0.033
betaExpansion = list()
meanInfectionExpansion = list()
with open(mathematicaFile, newline='') as csvfile:
    data = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
    for row in data:
        beta = row[0]
        V = row[1]
        betaExpansion.append(beta)
        if isIndependent:
            meanInfectionExpansion.append(calculateMeanInfectedIndependentFromV(V, gamma, beta, alpha, degreeHist, meanSimplexDegree))
        else:
            meanInfectionExpansion.append(calculateMeanInfectedDependent(V, gamma, beta, alpha, degreeHist))

plt.figure()
plt.plot(betaExpansion,meanInfectionExpansion)

betaTheory = np.linspace(min(betaExpansion), max(betaExpansion)+0.01, 45)

for betaVal in betaTheory:
    roots = simplexTheory.solveEquilibrium(gamma, betaVal, alpha, minDegree, maxDegree, meanSimplexDegree, degreeSequence=None, isIndependent=isIndependent, type=type, r=r, digits=digits)
    for root in roots:
        plt.scatter(betaVal, root, s=10, color='black')

plt.xlabel(r'$\beta$')
plt.ylabel('infected average')
plt.show()
