from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.polynomial import polyroots
import pickle

filename = 'equilibriaData10252019-225011'
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
betaList = np.linspace(0, 0.1, 1000)
print(betaCrit)
initialGuesses = np.linspace(0, 1, 2)
digitsOfPrecision = 5

def solveEquilbrium(gamma, betaList, alpha, minK, maxK, digits=4):
    equilibrium = list()
    for beta in betaList:
        roots = list()
        for initialGuess in initialGuesses:
            root, data, ier, msg = fsolve(uniformEquilibriumFunction, initialGuess,  args=(gamma, beta, alpha, minK, maxK), full_output=True)
            if ier == 1 and round(np.asscalar(root), digits) not in set(roots):

                roots.append(round(np.asscalar(root), digits))
        equilibrium.append(roots)
    return equilibrium

def uniformEquilibriumFunction(V, gamma, beta, alpha, minK, maxK):
    avgK = 0.5*(minK+maxK)
    frac = 1/(maxK-minK+1)
    sum = 0
    for k in range(minK, maxK+1):
        sum = sum + k**2*(beta + alpha*V)/(gamma + beta*k*V + alpha*k*V**2)
    return frac/avgK*sum - 1

def calculateAvgInfected(V, gamma, beta, alpha, minK, maxK):
    frac = 1/(maxK-minK+1)
    sum = 0
    for k in range(minK, maxK+1):
        sum = sum + k*(beta*V + alpha*V**2)/(gamma + beta*k*V + alpha*k*V**2)
    return frac*sum

equilibrium = solveEquilbrium(gamma, betaList, alpha, minK, maxK, digitsOfPrecision)

for i in range(len(betaList)):
    for root in equilibrium[i]:
        infectedAvg = calculateAvgInfected(root, gamma, betaList[i], alpha, minK, maxK)
        plt.scatter(betaList[i], infectedAvg, s=2, color='black')
plt.plot([np.min(betaList), np.max(betaList)], [0, 0], 'k--')
plt.ylim([-0.5, 1])
plt.xlim([np.min(betaList), np.max(betaList)])
plt.xlabel(r'$\beta$')
plt.ylabel('infected average')
plt.show()
