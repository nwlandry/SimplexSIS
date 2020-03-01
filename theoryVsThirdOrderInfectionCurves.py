import simplexTheory
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *
import csv

mathematicaFile = 'depData.csv'

gamma = 2
minDegree = 67
maxDegree = 450
isIndependent = False
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

betaTheory = np.linspace(min(betaExpansion)-0.001, max(betaExpansion)+0.01, 150)

for betaVal in betaTheory:
    roots = simplexTheory.solveEquilibrium(gamma, betaVal, alpha, minDegree, maxDegree, meanSimplexDegree, degreeSequence=None, isIndependent=isIndependent, type=type, r=r, digits=digits)
    for root in roots:
        plt.scatter(betaVal, root, s=5, color='black')

plt.xlabel(r'$\beta$')
plt.ylabel('infected average')
plt.show()
