import simplexTheory
import visualizeData
import simplexContagion
import simplexUtilities
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

plt.figure()
plt.hist(degreeSequence, bins=500)
plt.xlabel("degree")
plt.ylabel("Occurrences")
plt.show()

print(simplexUtilities.meanPowerOfDegree(degreeSequence, 4))
