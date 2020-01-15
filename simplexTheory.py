from scipy.optimize import fsolve, root
import matplotlib.pyplot as plt
import numpy as np

def solveEquilibrium(gamma, beta, alpha, minK, maxK, meanSimplexDegree, digits=4, isIndependent=False, type="power-law", r=4):
    if isIndependent:
        return solveIndependentEquilbrium(gamma, beta, alpha, minK, maxK, meanSimplexDegree, r=r, digits=digits, type=type)

    else:
        if type == "power-law":
            return solvePowerLawEquilbrium(gamma, beta, alpha, minK, maxK, r, digits)
        elif type == "uniform":
            return solveUniformEquilbrium(gamma, beta, alpha, minK, maxK, digits)
        else:
            print("invalid choice")

def solveUniformEquilbrium(gamma, beta, alpha, minK, maxK, digits=4):
    initialGuesses = np.linspace(0, 1, 5)
    roots = list()
    for initialGuess in initialGuesses:
        root, data, ier, msg = fsolve(uniformEquilibriumFunction, initialGuess,  args=(gamma, beta, alpha, minK, maxK), full_output=True)
        if ier == 1:
            avgInfectionPt = round(calculateAvgInfected(np.asscalar(root), gamma, beta, alpha, minK, maxK), digits)
            if avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
                roots.append(avgInfectionPt)
    return roots

def solvePowerLawEquilbrium(gamma, beta, alpha, minK, maxK, r, digits=4):
    initialGuesses = np.linspace(0, 1, 5)
    roots = list()
    for initialGuess in initialGuesses:
        root, data, ier, msg = fsolve(powerLawEquilibriumFunction, initialGuess,  args=(gamma, beta, alpha, minK, maxK, r), full_output=True)
        if ier == 1:
            avgInfectionPt = round(calculateAvgInfected(np.asscalar(root), gamma, beta, alpha, minK, maxK), digits)
            if avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
                roots.append(avgInfectionPt)
    return roots

def solveIndependentEquilbrium(gamma, beta, alpha, minK, maxK, meanSimplexDegree, r=4, digits=4, type="power-law"):
    initialGuesses = [[0, 0], [0.1, 0.1], [0.25, 0.25], [0.5, 0.5],[0.75, 0.75], [1, 1]]
    roots = list()
    for initialGuess in initialGuesses:
        result = root(independentEquilibriumFunction, initialGuess,  args=(gamma, beta, alpha, minK, maxK, meanSimplexDegree, r, type))
        avgInfectionPt = round(np.asscalar(result.x[0]), digits)
        if result.success and avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
            roots.append(avgInfectionPt)
    return roots

def uniformEquilibriumFunction(V, gamma, beta, alpha, minK, maxK):
    avgK = 0.5*(minK+maxK)
    frac = 1/(maxK-minK+1)
    sum = 0
    for k in range(minK, maxK+1):
        sum = sum + k**2*(beta + alpha*V)/(gamma + beta*k*V + alpha*k*V**2)
    return frac/avgK*sum - 1

def powerLawEquilibriumFunction(V, gamma, beta, alpha, minK, maxK, r):
    avgK = avgOfPowerLaw(minK, maxK, r)
    sum = 0
    for k in range(minK, maxK+1):
        sum = sum + truncatedPowerLaw(k, minK, maxK, r)*k**2*(beta + alpha*V)/(gamma + beta*k*V + alpha*k*V**2)
    return 1/avgK*sum - 1

def independentEquilibriumFunction(vars, gamma, beta, alpha, minK, maxK, meanSimplexDegree, r, networkDist="power-law"):
    # Add this calculation
    U = vars[0]
    V = vars[1]
    if networkDist == "power-law":
        avgK = avgOfPowerLaw(minK, maxK, r)
        sumV = 0
        sumU = 0
        for k in range(minK, maxK+1):
            sumU = sumU + truncatedPowerLaw(k, minK, maxK, r)*(k*beta*V + alpha*meanSimplexDegree*U**2)/(gamma + beta*k*V + alpha*meanSimplexDegree*U**2)
            sumV = sumV + truncatedPowerLaw(k, minK, maxK, r)*k*(k*beta*V + alpha*meanSimplexDegree*U**2)/(gamma + beta*k*V + alpha*meanSimplexDegree*U**2)

        return [sumU-U, 1/avgK*sumV - V]
    elif networkDist == "uniform":
        avgK = 0.5*(minK + maxK)
        frac = 1/(maxK-minK+1)
        sumV = 0
        sumU = 0
        for k in range(minK, maxK+1):
            sumU = sumU + (k*beta*V + alpha*meanSimplexDegree*U**2)/(gamma + beta*k*V + alpha*meanSimplexDegree*U**2)
            sumV = sumV + k*(k*beta*V + alpha*meanSimplexDegree*U**2)/(gamma + beta*k*V + alpha*meanSimplexDegree*U**2)

        return [frac*sumU-U, frac/avgK*sumV - V]
    else:
        print("invalid choice")
        return

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

def avgOfPowerLawEqn(minK, maxK, r, meanDeg):
    return (minK**(2-r)-maxK**(2-r))*(r-1)/((minK**(1-r)-maxK**(1-r))*(r-2)) - meanDeg

def meanPowerOfPowerLaw(minK, maxK, r, power):
    return (minK**(power+1-r)-maxK**(power+1-r))*(r-1)/((minK**(1-r)-maxK**(1-r))*(r-power-1))

def findHysteresisBoundaries(gamma, betaList, alphaList, minK, maxK, meanSimplexDegree, tolerance, isIndependent=False, type="uniform", r=4):
    # These depend on there being an odd number of points
    n = len(alphaList)
    m = len(betaList)
    leftBoundary = []
    rightBoundary = []
    alphaHys = []
    for i in range(n):
        left = max(betaList)
        right = 0
        isHysteresis = False
        for j in range(m):
            roots = solveEquilibrium(gamma, betaList[j], alphaList[i], minK, maxK, meanSimplexDegree, digits=4, isIndependent=isIndependent, type=type, r=r)
            if isIndependent:
                if len(roots) > 2:
                    isHysteresis = True
                    if betaList[j] < left:
                        left = betaList[j]

                    if betaList[j] > right:
                        right = betaList[j]
            else:
                if len(roots) > 1:
                    isHysteresis = True
                    if betaList[j] < left:
                        left = betaList[j]

                    if betaList[j] > right:
                        right = betaList[j]

        if isHysteresis:
            leftBoundary.append(left)
            rightBoundary.append(right)
            alphaHys.append(alphaList[i])

    # if type == "uniform":
    #     meanDegree = 0.5*(minK + maxK)
    #
    # elif type == "power-law":
    #     meanDegree = avgOfPowerLaw(minK, maxK, r)
    #
    # if not isIndependent:
    #     meanSimplexDegree = meanDegree

    return leftBoundary, rightBoundary, alphaHys
