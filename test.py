import simplexUtilities
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from datetime import *
from simplexTheory import *
import time


gamma = 2
minDegree = 1
maxDegree = 9
n = 10
isIndependent = True
type = "power-law"

r = 3
digits = 5

# degreeSequence = simplexUtilities.generatePowerLawDegreeSequence(n, minDegree, maxDegree, r)
# print(degreeSequence)

degreeSequence = [1,1,1,1,1,2,2,4,5]

# A = simplexUtilities.generateConfigModelAdjacency(degreeSequence)
# print(A)

simplexListDep, simplexIndices = simplexUtilities.generateConfigModelSimplexList(degreeSequence, 3)

print(simplexListDep)

meanSimplexDegree = sum(degreeSequence)/len(degreeSequence)
simplexListIndep, simplexIndices = simplexUtilities.generateUniformSimplexList(9, meanSimplexDegree, 3)

print(simplexListIndep)
