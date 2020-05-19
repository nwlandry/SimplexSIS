import simplexUtilities
import simplexTheory
import simplexVisualize
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from datetime import *

filename = 'equilibriaData05172020-072921'
#filename = 'Archive-Data/equilibriaData11112019-002636'
with open(filename, 'rb') as file:
    data = pickle.load(file)

gamma = data[0]
beta = data[1]
alpha = data[2]
equilibria = data[3]
degreeSequence = data[4]
isDegreeCorrelated = data[5]
type = data[6]
r = data[7]
print(type)
print(r)
print(isDegreeCorrelated)
print(simplexUtilities.meanPowerOfDegree(degreeSequence, 2))
print(min(degreeSequence))
