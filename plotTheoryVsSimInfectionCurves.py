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
filename = 'equilibriaData02092020-141544'
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
# alpha = data[0]
# beta = data[1]
# gamma = data[2]
# equilibria = data[5]

minDegree = 67
maxDegree = 450
isIndependent = True
type = "power-law"
meanSimplexDegree = 100
r = 4.001
digits = 5

index = 20
print(alpha[index])

#alphaCrit = meanPowerOfPowerLaw(minDegree, maxDegree, r, 3)/(meanPowerOfPowerLaw(minDegree, maxDegree, r, 1)**4)*gamma
#print(alphaCrit)
#visualizeData.plotTheoreticalInfectionCurves(gamma, beta, alpha[index], minDegree, maxDegree, meanSimplexDegree, digits=digits, isIndependent=isIndependent, type=type, r=r)

visualizeData.plotTheoreticalAndSimInfectionCurves(equilibria[index], gamma, beta, alpha[index], minDegree, maxDegree, meanSimplexDegree, degreeSequence=degreeSequence, isIndependent=isIndependent, type=type, r=r, digits=digits)
