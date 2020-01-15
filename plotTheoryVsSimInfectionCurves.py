import simplexTheory
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np


#filename = 'Poster/uniform_indep'
filename = 'Poster/power-law_r=4_indep'
#filename = 'Poster/power-law_r=4point5_indep'
#filename = 'Archive-Data/equilibriaData10252019-225011'
#filename = 'Archive-Data/equilibriaData11112019-002636'
with open(filename, 'rb') as file:
    data = pickle.load(file)

alpha = data[0]
beta = data[1]
gamma = data[2]
equilibria = data[5]

minK = 10
maxK = 200
isIndependent = True
type = "power-law"
meanSimplexDegree = 20
r = 3

index = 5
digits = 5
alphaCrit = 0.25
visualizeData.plotTheoreticalInfectionCurves(gamma, beta, alphaCrit, minK, maxK, meanSimplexDegree, digits=digits, isIndependent=isIndependent, type=type, r=r)

#visualizeData.plotTheoreticalAndSimInfectionCurves(equilibria[index], gamma, beta, alpha[index], minK, maxK, meanSimplexDegree, digits=digits, isIndependent=isIndependent, type=type, r=r)
