from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.polynomial import polyroots
import pickle
import simplexTheory

filename = 'equilibriaData11072019-102634'
with open(filename, 'rb') as file:
    data = pickle.load(file)
alpha = data[0]
beta = data[1]
equilibria = data[2]
plt.figure()
plt.plot(beta, equilibria[-1], 'o-', label=r"$\alpha=0.075$")


gamma = 2
alpha = 0.075
minK = 50
maxK = 100
avgK = 75
avgSquaredK = 5833.3
betaCrit = avgK/avgSquaredK*gamma
betaList = np.linspace(0, 0.1, 20)
print(betaCrit)
digitsOfPrecision = 5

# equilibriumV = simplexTheory.solveUniformEquilbrium(gamma, betaList, alpha, minK, maxK, digitsOfPrecision)
maxK = 1000
r = 4

equilibriumV = simplexTheory.solvePowerLawEquilbrium(gamma, betaList, alpha, minK, maxK, r, digitsOfPrecision)

for i in range(len(betaList)):
    for root in equilibriumV[i]:
        infectedAvg = simplexTheory.calculateAvgInfected(root, gamma, betaList[i], alpha, minK, maxK)
        plt.scatter(betaList[i], infectedAvg, s=2, color='black')
plt.plot([np.min(betaList), np.max(betaList)], [0, 0], 'k--')
plt.ylim([-0.5, 1])
plt.xlim([np.min(betaList), np.max(betaList)])
plt.xlabel(r'$\beta$')
plt.ylabel('infected average')
plt.show()
