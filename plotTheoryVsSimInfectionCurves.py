import simplexTheory
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *


#filename = 'Poster/uniform_indep'
filename = 'equilibriaData12262019-000110'
#filename = 'Poster/power-law_r=4_dep'
#filename = 'Archive-Data/equilibriaData10252019-225011'
#filename = 'Archive-Data/equilibriaData11112019-002636'
with open(filename, 'rb') as file:
    data = pickle.load(file)

alpha = data[0]
beta = data[1]
gamma = data[2]
equilibria = data[5]

minDegree = 13
maxDegree = 300
isIndependent = True
type = "power-law"
meanSimplexDegree = 20
r = 3.001

index = 6
digits = 5
alphaCrit = meanPowerOfPowerLaw(minDegree, maxDegree, r, 3)/(meanPowerOfPowerLaw(minDegree, maxDegree, r, 1)**4)*gamma
#print(alphaCrit)
visualizeData.plotTheoreticalInfectionCurves(gamma, beta, alphaCrit, minDegree, maxDegree, meanSimplexDegree, digits=digits, isIndependent=isIndependent, type=type, r=r)

#visualizeData.plotTheoreticalAndSimInfectionCurves(equilibria[index], gamma, beta, alpha[index], minDegree, maxDegree, meanSimplexDegree, digits=digits, isIndependent=isIndependent, type=type, r=r)
