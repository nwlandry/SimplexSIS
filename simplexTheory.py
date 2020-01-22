from scipy.optimize import fsolve, root
import matplotlib.pyplot as plt
import numpy as np

def solveEquilibrium(gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, digits=4, isIndependent=False, type="power-law", r=4):
    if isIndependent:
        # return solveIndependentEquilbrium(gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, r=r, digits=digits, type=type)
        return solveIndependentEquilbriumExpansion(gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, r=r, digits=digits, type=type)

    else:
        return solveDependentEquilbrium(gamma, beta, alpha, minDegree, maxDegree, r=r, digits=digits, type=type)

def solveDependentEquilbrium(gamma, beta, alpha, minDegree, maxDegree, r=4, digits=4, type="power-law"):
    initialGuesses = np.linspace(0, 1, 5)
    roots = list()
    for initialGuess in initialGuesses:
        root, data, ier, msg = fsolve(dependentEquilibriumFunction, initialGuess,  args=(gamma, beta, alpha, minDegree, maxDegree, r, type), full_output=True)
        if ier == 1:
            avgInfectionPt = round(calculateAvgInfected(np.asscalar(root), gamma, beta, alpha, minDegree, maxDegree), digits)
            if avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
                roots.append(avgInfectionPt)
    return roots

def solveIndependentEquilbriumExpansion(gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, r, digits=4, type="power-law"):
    initialGuesses = np.linspace(0, 1, 5)
    roots = list()
    for initialGuess in initialGuesses:
        root, data, ier, msg = fsolve(independentExpansionEquilibriumFunction, initialGuess,  args=(gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, r, type), full_output=True)
        if ier == 1:
            avgInfectionPt = round(calculateAvgInfected(np.asscalar(root), gamma, beta, alpha, minDegree, maxDegree), digits)
            if avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
                roots.append(avgInfectionPt)
    return roots

def solveIndependentEquilbrium(gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, r=4, digits=4, type="power-law"):
    initialGuesses = [[0, 0], [0.001, 0.001], [0.25, 0.25], [0.5, 0.5],[0.75, 0.75], [1, 1]]
    roots = list()
    for initialGuess in initialGuesses:
        result = root(independentEquilibriumFunction, initialGuess,  args=(gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, r, type))
        avgInfectionPt = round(np.asscalar(result.x[0]), digits)
        if result.success and avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
            roots.append(avgInfectionPt)
    return roots

def dependentEquilibriumFunction(V, gamma, beta, alpha, minDegree, maxDegree, r, networkDist):
    if networkDist == "power-law":
        meanDegree = avgOfPowerLaw(minDegree, maxDegree, r)
        sum = 0
        for k in range(minDegree, maxDegree+1):
            sum = sum + truncatedPowerLaw(k, minDegree, maxDegree, r)*k**2*(beta + alpha*V)/(gamma + beta*k*V + alpha*k*V**2)
        return 1/meanDegree*sum - 1

    elif networkDist == "uniform":
        meanDegree = 0.5*(minDegree+maxDegree)
        frac = 1/(maxDegree-minDegree+1)
        sum = 0
        for k in range(minDegree, maxDegree+1):
            sum = sum + k**2*(beta + alpha*V)/(gamma + beta*k*V + alpha*k*V**2)
        return frac/meanDegree*sum - 1

    else:
        print("invalid choice")
        return

def independentExpansionEquilibriumFunction(V, gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, r, networkDist):
    if networkDist == "power-law":
        meanDegree = avgOfPowerLaw(minDegree, maxDegree, r)
        meanSquaredDegree = meanPowerOfPowerLaw(minDegree, maxDegree, r, 2)
        c = meanDegree**4*meanSimplexDegree/meanSquaredDegree**2
        sum = 0
        for k in range(minDegree, maxDegree+1):
            sum = sum + truncatedPowerLaw(k, minDegree, maxDegree, r)*k*(k*beta + alpha*c*V)/(gamma + beta*k*V + alpha*c*V**2)
        return 1/meanDegree*sum - 1

    elif networkDist == "uniform":
        meanDegree = 0.5*(minDegree+maxDegree)
        meanSquaredDegree = 1./3*(maxDegree**2+minDegree*maxDegree+minDegree**2)
        frac = 1/(maxDegree-minDegree+1)
        c = meanDegree**4*meanSimplexDegree/meanSquaredDegree**2
        sum = 0
        for k in range(minDegree, maxDegree+1):
            sum = sum + k**2*(beta + alpha*c*V)/(gamma + beta*k*V + alpha*c*V**2)
        return frac/meanDegree*sum - 1

    else:
        print("invalid choice")
        return

def independentEquilibriumFunction(vars, gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, r, networkDist):
    # Add this calculation
    U = vars[0]
    V = vars[1]
    if networkDist == "power-law":
        meanDegree = avgOfPowerLaw(minDegree, maxDegree, r)
        sumV = 0
        sumU = 0
        for k in range(minDegree, maxDegree+1):
            sumU = sumU + truncatedPowerLaw(k, minDegree, maxDegree, r)*(k*beta*V + alpha*meanSimplexDegree*U**2)/(gamma + beta*k*V + alpha*meanSimplexDegree*U**2)
            sumV = sumV + truncatedPowerLaw(k, minDegree, maxDegree, r)*k*(k*beta*V + alpha*meanSimplexDegree*U**2)/(gamma + beta*k*V + alpha*meanSimplexDegree*U**2)

        return [sumU-U, 1/meanDegree*sumV - V]
    elif networkDist == "uniform":
        meanDegree = 0.5*(minDegree + maxDegree)
        frac = 1/(maxDegree-minDegree+1)
        sumV = 0
        sumU = 0
        for k in range(minDegree, maxDegree+1):
            sumU = sumU + (k*beta*V + alpha*meanSimplexDegree*U**2)/(gamma + beta*k*V + alpha*meanSimplexDegree*U**2)
            sumV = sumV + k*(k*beta*V + alpha*meanSimplexDegree*U**2)/(gamma + beta*k*V + alpha*meanSimplexDegree*U**2)

        return [frac*sumU-U, frac/meanDegree*sumV - V]
    else:
        print("invalid choice")
        return

def calculateAvgInfected(V, gamma, beta, alpha, minDegree, maxDegree):
    frac = 1/(maxDegree-minDegree+1)
    sum = 0
    for k in range(minDegree, maxDegree+1):
        sum = sum + k*(beta*V + alpha*V**2)/(gamma + beta*k*V + alpha*k*V**2)
    return frac*sum

def truncatedPowerLaw(k, minDegree, maxDegree, r):
    return (r-1)/(minDegree**(1-r)-maxDegree**(1-r))*k**(-r)

def avgOfPowerLaw(minDegree, maxDegree, r):
    return (minDegree**(2-r)-maxDegree**(2-r))*(r-1)/((minDegree**(1-r)-maxDegree**(1-r))*(r-2))

def avgOfPowerLawEqn(minDegree, maxDegree, r, meanDeg):
    return (minDegree**(2-r)-maxDegree**(2-r))*(r-1)/((minDegree**(1-r)-maxDegree**(1-r))*(r-2)) - meanDeg

def meanPowerOfPowerLaw(minDegree, maxDegree, r, power):
    return (minDegree**(power+1-r)-maxDegree**(power+1-r))*(r-1)/((minDegree**(1-r)-maxDegree**(1-r))*(r-power-1))

def findHysteresisBoundaries(gamma, betaList, alphaList, minDegree, maxDegree, meanSimplexDegree, tolerance, isIndependent=False, type="uniform", r=4):
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
            roots = solveEquilibrium(gamma, betaList[j], alphaList[i], minDegree, maxDegree, meanSimplexDegree, digits=4, isIndependent=isIndependent, type=type, r=r)
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
    #     meanDegree = 0.5*(minDegree + maxDegree)
    #
    # elif type == "power-law":
    #     meanDegree = avgOfPowerLaw(minDegree, maxDegree, r)
    #
    # if not isIndependent:
    #     meanSimplexDegree = meanDegree

    return leftBoundary, rightBoundary, alphaHys

### Additional functions
def solveUniformEquilbrium(gamma, beta, alpha, minDegree, maxDegree, digits=4):
    initialGuesses = np.linspace(0, 1, 5)
    roots = list()
    for initialGuess in initialGuesses:
        root, data, ier, msg = fsolve(uniformEquilibriumFunction, initialGuess,  args=(gamma, beta, alpha, minDegree, maxDegree), full_output=True)
        if ier == 1:
            avgInfectionPt = round(calculateAvgInfected(np.asscalar(root), gamma, beta, alpha, minDegree, maxDegree), digits)
            if avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
                roots.append(avgInfectionPt)
    return roots

def solvePowerLawEquilbrium(gamma, beta, alpha, minDegree, maxDegree, r, digits=4):
    initialGuesses = np.linspace(0, 1, 5)
    roots = list()
    for initialGuess in initialGuesses:
        root, data, ier, msg = fsolve(powerLawEquilibriumFunction, initialGuess,  args=(gamma, beta, alpha, minDegree, maxDegree, r), full_output=True)
        if ier == 1:
            avgInfectionPt = round(calculateAvgInfected(np.asscalar(root), gamma, beta, alpha, minDegree, maxDegree), digits)
            if avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
                roots.append(avgInfectionPt)
    return roots


def uniformEquilibriumFunction(V, gamma, beta, alpha, minDegree, maxDegree):
    meanDegree = 0.5*(minDegree+maxDegree)
    frac = 1/(maxDegree-minDegree+1)
    sum = 0
    for k in range(minDegree, maxDegree+1):
        sum = sum + k**2*(beta + alpha*V)/(gamma + beta*k*V + alpha*k*V**2)
    return frac/meanDegree*sum - 1

def powerLawEquilibriumFunction(V, gamma, beta, alpha, minDegree, maxDegree, r):
    meanDegree = avgOfPowerLaw(minDegree, maxDegree, r)
    sum = 0
    for k in range(minDegree, maxDegree+1):
        sum = sum + truncatedPowerLaw(k, minDegree, maxDegree, r)*k**2*(beta + alpha*V)/(gamma + beta*k*V + alpha*k*V**2)
    return 1/meanDegree*sum - 1
