import simplexTheory
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *

# gamma = 2
# isIndependent = False
# type = "power-law"
# minDegree = 50
# maxDegreeList = np.linspace(100, 1000, 19)
# rList = np.linspace(2.5,4.0,16)
#
# digits = 4
# tolerance = 0.001
# numBetaPoints = 30
# minAlpha = 0.02
# maxAlpha = 0.1
# m = np.size(maxDegreeList,0)
# n = np.size(rList,0)
# alphaCritGrid = np.zeros([m, n])
#
# for i in range(m):
#     for j in range(n):
#         degreeHist = generateTheoreticalDegreeHist(minDegree, int(maxDegreeList[i]), type, r=rList[j])
#
#         meanDegree = sum([k*prob for k, prob in degreeHist])
#         meanSquaredDegree = sum([k**2*prob for k, prob in degreeHist])
#         meanSimplexDegree = meanDegree
#
#         betaTheory = np.linspace(0.5*meanDegree/meanSquaredDegree*gamma, 1.5*meanDegree/meanSquaredDegree*gamma, 9)
#
#         alphaCrit = calculateTheoreticalCriticalAlpha(gamma, betaTheory, minAlpha, maxAlpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, option="infinity", digits=4, tolerance=tolerance)
#         alphaCritGrid[i][j] = alphaCrit
#
#         print("kmax="+str(maxDegreeList[i]) + ", r=" + str(rList[j]))
#
# xMin = np.min(rList)
# xMax = np.max(rList)
# yMin = np.min(maxDegreeList)
# yMax = np.max(maxDegreeList)

filename = "alphaCrit03292020-125621"
with open(filename, 'rb') as file:
    data = pickle.load(file)

xMin = data[0]
xMax = data[1]
yMin = data[2]
yMax = data[3]
alphaCritGrid = data[4]

plt.figure()
c = plt.imshow(np.flipud(alphaCritGrid), interpolation="none", cmap="Reds", extent=[xMin, xMax, yMin, yMax], aspect="auto")
cbar = plt.colorbar(c)
cbar.set_label(r"$\alpha_{crit}$", rotation=90)
plt.xlabel("Power-Law Exponent")
plt.ylabel("Maximum degree")
plt.plot()
plt.show()
