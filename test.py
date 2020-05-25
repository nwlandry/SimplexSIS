import simplexTheory
import simplexVisualize
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *
import multiprocessing as mp
from datetime import datetime
import time
import os
import math
import scipy.io
from scipy.optimize import fsolve

filename = "Paper Figures/alphaCritUncorrelated"
#filename = "alphaCrit05122020-004430"
with open(filename, 'rb') as file:
    data = pickle.load(file)
#
# def f(beta3, k1, k2, k3, k4, k5, gamma):
#     return 4*k1**11*beta3**3 - 4*k1**6*k2**2*beta3**2*gamma - 16*k1**7*k3*beta3**2*gamma + 8*k1**5*k2*k3*beta3**2*gamma - 24*k1*k2*k3**2*beta3*gamma**2 + 4*k1**2*k3*(5*k2**2+k4)*beta3*gamma**2 + 4*k1**3*(3*k3**2 - 2*k2*k4)*beta3*gamma**2 + 4*k1**4*k5*beta3*gamma**2 + gamma**2*(4*k3**3*beta3 + k4**2*gamma - 4*k3*k5*gamma)


def criticalBeta3(degreeHist, isDegreeCorrelated, gamma):

    k1 = computeMeanPowerOfDegreeFromHist(degreeHist, 1)
    k2 = computeMeanPowerOfDegreeFromHist(degreeHist, 2)
    k3 = computeMeanPowerOfDegreeFromHist(degreeHist, 3)
    k4 = computeMeanPowerOfDegreeFromHist(degreeHist, 4)
    k5 = computeMeanPowerOfDegreeFromHist(degreeHist, 5)

    if isDegreeCorrelated:
        return k3*k1**2/k2**3*gamma

    else:

        def a2(beta3, k1, k2, k3, gamma):
            return (beta3*k1**5 - gamma*k1*k3)/(k2**2*gamma)

        def a3(beta3, k1, k2, k3, k4, gamma):
            return k1**2*(gamma*k4 - 2*beta3*k1*(2*k1**2*k2 - k1*k3))/(k2**3*gamma)

        def a4(beta3, k1, k2, k3, k4, k5, gamma):
            return -k1**2*(gamma**2*k1*k5 - gamma*beta3*(3*k1**4*k3 + 4*k1**2*k2*(k1*k2 - k3) + k1*(k3 - k1*k2)**2) + k1**8*beta3**2)/(k2**4*gamma**2)

        def f1(beta3, k1, k2, k3, k4, k5, gamma):
            return a3(beta3, k1, k2, k3, k4, gamma)**2 - 4*a2(beta3, k1, k2, k3, gamma)*a4(beta3, k1, k2, k3, k4, k5, gamma)
        def f2(beta3, k1, k2, k3, k4, k5, gamma):
            return a3(beta3, k1, k2, k3, k4, gamma)**2 - a3(beta3, k1, k2, k3, k4, gamma)**2 - 4*a2(beta3, k1, k2, k3, gamma)*a4(beta3, k1, k2, k3, k4, k5, gamma)

        b1 = fsolve(f1, 0.1, args=(k1, k2, k3, k4, k5, gamma))

        # saddle-node bifurcation
        if a4(b1, k1, k2, k3, k4, k5, gamma) < 0 and -a3(b1, k1, k2, k3, k4, gamma)/(2*a4(b1, k1, k2, k3, k4, k5, gamma)) > 0 and -a3(b1, k1, k2, k3, k4, gamma)/(2*a4(b1, k1, k2, k3, k4, k5, gamma)) < 1:
            return b1
        else:
            return k3/k1**4*gamma

xMin = data[0]
xMax = data[1]
yMin = data[2]
yMax = data[3]
alphaCritGrid = data[4]

gamma = 2
isDegreeCorrelated = False
type = "power-law"
minDegree = 50
maxDegreeList = np.linspace(100, 1000, 37)
rList = np.linspace(2.5,4.0,31)

digits = 4
m = np.size(maxDegreeList,0)
n = np.size(rList,0)

firstOrderAlphaCritGrid = np.zeros([m,n])
higherOrderAlphaCritGrid = np.zeros([m,n])

for i in range(m):
    for j in range(n):
        degreeHist = generateTheoreticalDegreeHist(minDegree, int(maxDegreeList[i]), type, r=rList[j])
        higherOrderAlphaCritGrid[i,j] = criticalBeta3(degreeHist, isDegreeCorrelated, gamma)


# plt.figure()
# plt.subplot(121)
# x = np.linspace(xMin, xMax, np.size(alphaCritGrid, axis=1))
# y = np.linspace(yMin, yMax, np.size(alphaCritGrid, axis=0))
# c = plt.contour(x, y, higherOrderAlphaCritGrid, 10)
# plt.clabel(c, inline=1, fontsize=10)
#
# plt.subplot(122)
# c = plt.contour(x, y, alphaCritGrid, 10)
# plt.clabel(c, inline=1, fontsize=10)
# plt.show()
#
# plt.figure()
# x = np.linspace(xMin, xMax, np.size(alphaCritGrid, axis=1))
# y = np.linspace(yMin, yMax, np.size(alphaCritGrid, axis=0))
# c = plt.contour(x, y, np.divide(higherOrderAlphaCritGrid-alphaCritGrid,alphaCritGrid), 15)
# plt.clabel(c, inline=1, fontsize=10)
# plt.xlabel("Power-Law Exponent")
# plt.ylabel("Maximum degree")
# plt.show()

#
# plt.figure()
# c = plt.imshow(np.flipud(higherOrderAlphaCritGrid), interpolation="None", cmap="Reds", extent=[xMin, xMax, yMin, yMax], aspect="auto")
# cbar = plt.colorbar(c)
# cbar.set_label(r"$\alpha_{crit}$", rotation=90)
# plt.xlabel("Power-Law Exponent")
# plt.ylabel("Maximum degree")
# plt.show()
#

plt.figure()
c = plt.imshow(np.flipud(np.divide(higherOrderAlphaCritGrid-alphaCritGrid,alphaCritGrid)), interpolation="spline16", cmap="Reds", extent=[xMin, xMax, yMin, yMax], aspect="auto")
cbar = plt.colorbar(c)
cbar.set_label("Relative error", rotation=90)
plt.xlabel("Power-Law Exponent")
plt.ylabel("Maximum degree")
plt.show()
