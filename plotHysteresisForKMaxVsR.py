import simplexTheory
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *

gamma = 2
beta = 0.02
alpha = 0.05
isIndependent = True
type = "power-law"
minDegree = 50
maxDegreeList = np.linspace(100, 500, 9)
rList = np.linspace(3.0,4.0,10)

digits = 4
m = np.size(maxDegreeList,0)
n = np.size(rList,0)
hysteresis = np.zeros([m, n])

for i in range(m):
    for j in range(n):
        degreeHist = generateTheoreticalDegreeHist(minDegree, int(maxDegreeList[i]), type, r=rList[j])
        meanSimplexDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 1)
        roots = solveEquilibrium(gamma, beta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, digits=digits)
        if len(roots) == 3:
            hysteresis[i][j] = max(roots)-min(roots)


xMin = np.min(rList)
xMax = np.max(rList)
yMin = np.min(maxDegreeList)
yMax = np.max(maxDegreeList)


plt.figure()
xMin = min(rList)
c = plt.imshow(np.flipud(hysteresis), interpolation="spline16", cmap="Reds", extent=[xMin, xMax, yMin, yMax], aspect="auto")
cbar = plt.colorbar(c)
cbar.set_label('Distance between fixed point solutions', rotation=90)
plt.xlabel("Power-Law Exponent")
plt.ylabel("Maximum degree")
plt.title(r"Hysteresis for $\gamma=$"+str(gamma)+r", $\beta=$"+str(beta)+r", $\alpha=$"+str(alpha))
plt.plot()
plt.show()
