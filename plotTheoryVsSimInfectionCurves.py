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
filename = 'equilibriaData02292020-015156'
#filename = 'Archive-Data/equilibriaData10252019-225011'
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

degreeHist = simplexTheory.degreeSequenceToHist(degreeSequence)

index = 0
print(alpha[index])

visualizeData.plotTheoreticalAndSimInfectionCurves(equilibria[index], gamma, beta, alpha[index], degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, digits=digits)
