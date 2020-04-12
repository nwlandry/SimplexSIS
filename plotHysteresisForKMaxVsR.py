import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *

gamma = 1.0/14.0
R0 = 2.6
fracIndivid = 0.5

isIndependent = True
type = "power-law"
minDegree = 50
maxDegreeList = np.linspace(100, 500, 4)
rList = np.linspace(3.0,4.0,4)

digits = 4
m = np.size(maxDegreeList,0)
n = np.size(rList,0)
hysteresis = np.zeros([m, n])

for i in range(m):
    for j in range(n):
        degreeHist = generateTheoreticalDegreeHist(minDegree, int(maxDegreeList[i]), type, r=rList[j])
        meanDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 1)
        meanSquaredDegree = computeMeanPowerOfDegreeFromHist(degreeHist, 2)
        meanSimplexDegree = meanDegree

        beta = R0*meanDegree/meanSquaredDegree*gamma
        alpha = 2*R0*meanDegree/meanSquaredDegree*gamma
        print(alpha)
        # minAlpha = 0
        # maxAlpha = 0.06
        # minBeta = 0*meanDegree/meanSquaredDegree*gamma
        # maxBeta = 1.5*meanDegree/meanSquaredDegree*gamma
        # numBetaPoints = 75
        # tolerance = 0.001

        # alphaCrit = calculateTheoreticalCriticalAlpha(gamma, minBeta, maxBeta, minAlpha, maxAlpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, option="infinity", digits=digits, tolerance=tolerance)
        # print(alphaCrit)
        roots = solveEquilibrium(gamma, beta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, digits=digits)
        if len(roots) == 3:
            hysteresis[i][j] = max(roots)


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
plt.title(r"Hysteresis for $\gamma=$"+str(round(gamma,3))+r", $\beta=$"+str(round(beta,3))+r", $\alpha=$"+str(round(alpha,3)))
plt.plot()
plt.show()
