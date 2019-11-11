from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import numpy as np

def solveUniformEquilbrium(gamma, betaList, alpha, minK, maxK, digits=4):
    initialGuesses = np.linspace(0, 1, 2)
    equilibriumV = list()
    for beta in betaList:
        roots = list()
        for initialGuess in initialGuesses:
            root, data, ier, msg = fsolve(uniformEquilibriumFunction, initialGuess,  args=(gamma, beta, alpha, minK, maxK), full_output=True)
            if ier == 1 and round(np.asscalar(root), digits) not in set(roots):

                roots.append(round(np.asscalar(root), digits))
        equilibriumV.append(roots)
    return equilibriumV

def solvePowerLawEquilbrium(gamma, betaList, alpha, minK, maxK, r, digits=4):
    initialGuesses = np.linspace(0, 1, 2)
    equilibriumV = list()
    for beta in betaList:
        roots = list()
        for initialGuess in initialGuesses:
            root, data, ier, msg = fsolve(powerLawEquilibriumFunction, initialGuess,  args=(gamma, beta, alpha, minK, maxK, r), full_output=True)
            if ier == 1 and round(np.asscalar(root), digits) not in set(roots):
                roots.append(round(np.asscalar(root), digits))
        equilibriumV.append(roots)
    return equilibriumV

def uniformEquilibriumFunction(V, gamma, beta, alpha, minK, maxK):
    avgK = 0.5*(minK+maxK)
    frac = 1/(maxK-minK+1)
    sum = 0
    for k in range(minK, maxK+1):
        sum = sum + k**2*(beta + alpha*V)/(gamma + beta*k*V + alpha*k*V**2)
    return frac/avgK*sum - 1

def powerLawEquilibriumFunction(V, gamma, beta, alpha, minK, maxK, r):
    # Add this calculation
    avgK = avgOfPowerLaw(minK, maxK, r)
    sum = 0
    for k in range(minK, maxK+1):
        sum = sum + truncatedPowerLaw(k, minK, maxK, r)*k**2*(beta + alpha*V)/(gamma + beta*k*V + alpha*k*V**2)
    return 1/avgK*sum - 1

def calculateAvgInfected(V, gamma, beta, alpha, minK, maxK):
    frac = 1/(maxK-minK+1)
    sum = 0
    for k in range(minK, maxK+1):
        sum = sum + k*(beta*V + alpha*V**2)/(gamma + beta*k*V + alpha*k*V**2)
    return frac*sum

def truncatedPowerLaw(k, minK, maxK, r):
    return (r-1)/(minK**(1-r)-maxK**(1-r))*k**(-r)

def avgOfPowerLaw(minK, maxK, r):
    return (minK**(2-r)-maxK**(2-r))*(r-1)/((minK**(1-r)-maxK**(1-r))*(r-2))
