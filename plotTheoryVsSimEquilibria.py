import pickle
import matplotlib.pyplot as plt
import numpy as np
import simplexUtilities
import simplexVisualize
import simplexTheory

filename = 'equilibriaData06062020-003754'
with open(filename, 'rb') as file:
    data = pickle.load(file)
gamma = data[0]
beta = data[1]
alpha = data[2]
equilibria = data[3]
degreeSequence = data[4]
isDegreeCorrelated = data[5]
degreeDistType = data[6]
exponent = data[7]
meanSimplexDegree = data[8]
print(data[6])
print(data[7])
print(data[5])
index = 0
print(alpha[index])

degreeHist = simplexTheory.degreeSequenceToHist(degreeSequence)

simplexVisualize.plotTheoreticalAndSimInfectionCurves(equilibria[index], gamma, beta, alpha[index], degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, numTheoryPoints=100, digits=5)
