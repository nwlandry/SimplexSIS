import simplexTheory
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *
import csv

mathematicaFile = 'indepData.csv'

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

alpha = 0.023
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

betaTheory = np.linspace(min(betaExpansion), max(betaExpansion), 200)

for betaVal in betaTheory:
    roots = simplexTheory.solveEquilibrium(gamma, betaVal, alpha, minDegree, maxDegree, meanSimplexDegree, degreeSequence=None, isIndependent=isIndependent, type=type, r=r, digits=digits)
    for root in roots:
        plt.scatter(betaVal, root, s=5, color='black')

meanDegree = sum([k*prob for k, prob in degreeHist])
meanSquaredDegree = sum([k**2*prob for k, prob in degreeHist])
print(meanDegree)
print(meanSquaredDegree)

plt.scatter(98.5976/11499.9*2, 0, s=50, color='red')
plt.scatter(meanDegree/meanSquaredDegree*2, 0, s=50, color='green')
plt.xlabel(r'$\beta$')
plt.ylabel('infected average')
plt.ylim([0, 0.06])
plt.show()
