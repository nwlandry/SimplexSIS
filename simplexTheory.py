from scipy.optimize import fsolve, root
import matplotlib.pyplot as plt
import numpy as np

def getPhase(gamma, beta, alpha, degreeHist, meanSimplexDegree=None, isIndependent=False, majorityVote=True, healing=False, digits=4):
    return len(solveEquilibrium(gamma, beta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, majorityVote=majorityVote, healing=healing, digits=digits))

def solveEquilibrium(gamma, beta, alpha, degreeHist, meanSimplexDegree=None, isIndependent=False, majorityVote=True, healing=False, digits=4):
    if meanSimplexDegree is None:
        meanSimplexDegree = generateMeanDegreeFromHist(degreeHist)

    if majorityVote:
        if isIndependent:
            return solveIndependentEquilbriumMajorityVote(gamma, beta, alpha, degreeHist, meanSimplexDegree, healing=healing, digits=digits)
        else:
            return solveDependentEquilbriumMajorityVote(gamma, beta, alpha, degreeHist, healing=healing, digits=digits)
    else:
        if isIndependent:
            return solveIndependentEquilbriumIndividual(gamma, beta, alpha, degreeHist, meanSimplexDegree, healing=healing, digits=digits)
        else:
            return solveDependentEquilbriumIndividual(gamma, beta, alpha, degreeHist, healing=healing, digits=digits)

def solveDependentEquilbriumMajorityVote(gamma, beta, alpha, degreeHist, healing=False, digits=4):
    initialGuesses = np.linspace(0, 0.6, 10)
    roots = list()
    for initialGuess in initialGuesses:
        root, data, ier, msg = fsolve(dependentEquilibriumFunctionMajorityVote, initialGuess,  args=(gamma, beta, alpha, degreeHist, healing), full_output=True)
        if ier == 1:
            avgInfectionPt = round(calculateMeanInfectedDependentMajorityVote(np.asscalar(root), gamma, beta, alpha, degreeHist, healing), digits)
            if avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
                roots.append(avgInfectionPt)
    return roots

def solveIndependentEquilbriumMajorityVote(gamma, beta, alpha, degreeHist, meanSimplexDegree, healing=False, digits=4):
    initialGuesses = [[val,val] for val in np.linspace(0, 0.6, 10)]
    roots = list()
    for initialGuess in initialGuesses:
        result = root(independentEquilibriumFunctionMajorityVote, initialGuess,  args=(gamma, beta, alpha, degreeHist, meanSimplexDegree, healing))
        avgInfectionPt = round(np.asscalar(result.x[0]), digits)
        if result.success and avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
            roots.append(avgInfectionPt)
    return roots

def solveDependentEquilbriumIndividual(gamma, beta, alpha, degreeHist, healing=False, digits=4):
    initialGuesses = np.linspace(0, 0.6, 10)
    roots = list()
    for initialGuess in initialGuesses:
        root, data, ier, msg = fsolve(dependentEquilibriumFunctionIndividual, initialGuess,  args=(gamma, beta, alpha, degreeHist, healing), full_output=True)
        if ier == 1:
            avgInfectionPt = round(calculateMeanInfectedDependentIndividual(np.asscalar(root), gamma, beta, alpha, degreeHist, healing), digits)
            if avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
                roots.append(avgInfectionPt)
    return roots

def solveIndependentEquilbriumIndividual(gamma, beta, alpha, degreeHist, meanSimplexDegree, healing=False, digits=4):
    initialGuesses = [[val,val] for val in np.linspace(0, 0.6, 10)]
    roots = list()
    for initialGuess in initialGuesses:
        result = root(independentEquilibriumFunctionIndividual, initialGuess,  args=(gamma, beta, alpha, degreeHist, meanSimplexDegree, healing))
        avgInfectionPt = round(np.asscalar(result.x[0]), digits)
        if result.success and avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
            roots.append(avgInfectionPt)
    return roots

def dependentEquilibriumFunctionMajorityVote(V, gamma, beta, alpha, degreeHist, healing=False):
    meanDegree = 0
    sum = 0
    if healing:
        sign = -1
    else:
        sign = 1
    for degreeInfo in degreeHist:
        degree = degreeInfo[0]
        prob = degreeInfo[1]
        sum = sum + prob*degree**2*(beta*V + sign*alpha*V**2)/(gamma + beta*degree*V + sign*alpha*degree*V**2)
        meanDegree = meanDegree + prob*degree
    return 1/meanDegree*sum - V

def independentEquilibriumFunctionMajorityVote(vars, gamma, beta, alpha, degreeHist, meanSimplexDegree, healing=False):
    # Add this calculation
    U = vars[0]
    V = vars[1]
    meanDegree = 0
    sumU = 0
    sumV = 0
    if healing:
        sign = -1
    else:
        sign = 1
    for degreeInfo in degreeHist:
        degree = degreeInfo[0]
        prob = degreeInfo[1]
        sumU = sumU + prob*(degree*beta*V + sign*alpha*meanSimplexDegree*U**2)/(gamma + beta*degree*V + sign*alpha*meanSimplexDegree*U**2)
        sumV = sumV + prob*degree*(degree*beta*V + sign*alpha*meanSimplexDegree*U**2)/(gamma + beta*degree*V + sign*alpha*meanSimplexDegree*U**2)
        meanDegree = meanDegree + prob*degree

    return [sumU-U, 1/meanDegree*sumV - V]

def dependentEquilibriumFunctionIndividual(V, gamma, beta, alpha, degreeHist, healing=False):
    meanDegree = 0
    sum = 0
    if healing:
        sign = -1
    else:
        sign = 1
    for degreeInfo in degreeHist:
        degree = degreeInfo[0]
        prob = degreeInfo[1]
        sum = sum + prob*degree**2*((beta + sign*2*alpha)*V - sign*alpha*V**2)/(gamma + (beta + sign*2*alpha)*degree*V - sign*alpha*degree*V**2)
        meanDegree = meanDegree + prob*degree
    return 1/meanDegree*sum - V

def independentEquilibriumFunctionIndividual(vars, gamma, beta, alpha, degreeHist, meanSimplexDegree, healing=False):
    # Add this calculation
    U = vars[0]
    V = vars[1]
    meanDegree = 0
    sumU = 0
    sumV = 0
    if healing:
        sign = -1
    else:
        sign = 1
    for degreeInfo in degreeHist:
        degree = degreeInfo[0]
        prob = degreeInfo[1]
        sumU = sumU + prob*(degree*beta*V + sign*2*alpha*meanSimplexDegree*U - sign*alpha*meanSimplexDegree*U**2)/(gamma + beta*degree*V + sign*2*alpha*meanSimplexDegree*U - sign*alpha*meanSimplexDegree*U**2)
        sumV = sumV + prob*degree*(degree*beta*V + sign*2*alpha*meanSimplexDegree*U - sign*alpha*meanSimplexDegree*U**2)/(gamma + beta*degree*V + sign*2*alpha*meanSimplexDegree*U - sign*alpha*meanSimplexDegree*U**2)
        meanDegree = meanDegree + prob*degree

    return [sumU-U, 1/meanDegree*sumV - V]

def degreeSequenceToHist(degreeSequence):
    degreeHist = list()
    sortedDegreeSequence = sorted(degreeSequence)
    n = len(degreeSequence)
    currentDegree = 0
    for degree in sortedDegreeSequence:
        if degree != currentDegree:
            currentDegree = degree
            degreeHist.append([degree, 1])
        else:
            degreeHist[-1][1] = degreeHist[-1][1] + 1
    for i in range(len(degreeHist)):
        degreeHist[i][1] = degreeHist[i][1]/n
    return degreeHist

def generateTheoreticalDegreeHist(minDegree, maxDegree, networkDist, r=4):
    degreeHist = list()
    for degree in range(minDegree, maxDegree+1):
        if networkDist == "power-law":
            prob = truncatedPowerLaw(degree, minDegree, maxDegree, r)
        elif networkDist == "uniform":
            prob = 1.0/(maxDegree-minDegree+1)
        else:
            print("Not a valid option")
        degreeHist.append([degree, prob])
    totalProb = sum([prob for k, prob in degreeHist])
    if totalProb != 1:
        for item in degreeHist:
            item[1] = item[1]/totalProb
    return degreeHist

def computeMeanPowerOfDegreeFromHist(hist, power):
    meanDegree = 0
    for degree, probability in hist:
        meanDegree = meanDegree + degree**power*probability
    return meanDegree

def calculateMeanInfectedDependentMajorityVote(V, gamma, beta, alpha, degreeHist, healing=False):
    meanInfected = 0
    if healing:
        sign = -1
    else:
        sign = 1
    for degreeInfo in degreeHist:
        probability = degreeInfo[1]
        degree = degreeInfo[0]
        meanInfected = meanInfected + probability*degree*(beta*V + sign*alpha*V**2)/(gamma + beta*degree*V + sign*alpha*degree*V**2)
    return meanInfected

def calculateMeanInfectedDependentIndividual(V, gamma, beta, alpha, degreeHist, healing=False):
    meanInfected = 0
    if healing:
        sign = -1
    else:
        sign = 1
    for degreeInfo in degreeHist:
        probability = degreeInfo[1]
        degree = degreeInfo[0]
        meanInfected = meanInfected + probability*degree*((beta + sign*2*alpha)*V - sign*alpha*V**2)/(gamma + (beta + sign*2*alpha)*degree*V - sign*alpha*degree*V**2)
    return meanInfected

def calculateMeanInfectedIndependentFromV(V, gamma, beta, alpha, degreeHist, meanSimplexDegree):
    initialGuess = 0.5
    meanInfected = fsolve(fV, initialGuess,  args=(V, gamma, beta, alpha, degreeHist, meanSimplexDegree))
    return meanInfected

def fV(U, V, gamma, beta, alpha, degreeHist, meanSimplexDegree):
    sumU = 0
    for degreeInfo in degreeHist:
        probability = degreeInfo[1]
        degree = degreeInfo[0]
        sumU = sumU + probability*(beta*degree*V + alpha*meanSimplexDegree*U**2)/(gamma + beta*degree*V + alpha*meanSimplexDegree*U**2)

    return sumU-U

def truncatedPowerLaw(k, minDegree, maxDegree, r):
    return (r-1)/(minDegree**(1-r)-maxDegree**(1-r))*k**(-r)

def avgOfPowerLaw(minDegree, maxDegree, r):
    return (minDegree**(2-r)-maxDegree**(2-r))*(r-1)/((minDegree**(1-r)-maxDegree**(1-r))*(r-2))

def avgOfPowerLawEqn(minDegree, maxDegree, r, meanDeg):
    return (minDegree**(2-r)-maxDegree**(2-r))*(r-1)/((minDegree**(1-r)-maxDegree**(1-r))*(r-2)) - meanDeg

def meanPowerOfPowerLaw(minDegree, maxDegree, r, power):
    return (minDegree**(power+1-r)-maxDegree**(power+1-r))*(r-1)/((minDegree**(1-r)-maxDegree**(1-r))*(r-power-1))

def calculateTheoreticalHysteresisVisually(gamma, minBeta, maxBeta, alpha, degreeHist, meanSimplexDegree=None, isIndependent=False, option="infinity", digits=4, tolerance=0.0001):
    if meanSimplexDegree == None:
        meanSimplexDegree = sum([k*prob for k, prob in degreeHist])

    plt.figure()
    minRoots = solveEquilibrium(gamma, minBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, digits=digits)

    maxRoots = solveEquilibrium(gamma, maxBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, digits=digits)

    plt.scatter(minBeta, len(minRoots))
    plt.scatter(maxBeta, len(maxRoots))
    hysteresis = 0
    if (max(minRoots) < tolerance and max(maxRoots) > tolerance) or len(minRoots) == 3:
        while maxBeta - minBeta > tolerance:
            newBeta = 0.5*(minBeta + maxBeta)
            newRoots = solveEquilibrium(gamma, newBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, digits=digits)
            plt.scatter(newBeta,len(newRoots))
            if len(newRoots) == 3:
                if option == "infinity":
                    if max(newRoots) >= hysteresis: # min of the roots is always 0
                        hysteresis = max(newRoots)
                        minBeta = newBeta # Because the epidemic fraction increases with increasing beta
            elif len(newRoots) == 2:
                testBeta = newBeta - tolerance/2.0
                testRoots = solveEquilibrium(gamma, testBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, digits=digits) # perturb the solution to see if this is at the end of the branch or not
                if len(testRoots) == 1:
                    minBeta = newBeta
                    # if max(newRoots) >= hysteresis and max(newRoots)>0.1: # min of the roots is always 0
                    #     hysteresis = max(newRoots)
                else:
                    maxBeta = newBeta
            else: # if just the zero solution
                minBeta = newBeta
        plt.show()
        return hysteresis
    else:
        return float("nan")

def calculateTheoreticalHysteresis(gamma, minBeta, maxBeta, alpha, degreeHist, meanSimplexDegree=None, isIndependent=False, option="infinity", digits=4, tolerance=0.0001):
    if meanSimplexDegree == None:
        meanSimplexDegree = sum([k*prob for k, prob in degreeHist])

    minRoots = solveEquilibrium(gamma, minBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, digits=digits)

    maxRoots = solveEquilibrium(gamma, maxBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, digits=digits)
    hysteresis = 0
    if (max(minRoots) < tolerance and max(maxRoots) > tolerance) or len(minRoots) == 3:
        while maxBeta - minBeta > tolerance:
            newBeta = 0.5*(minBeta + maxBeta)
            newRoots = solveEquilibrium(gamma, newBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, digits=digits)

            if len(newRoots) == 3:
                if option == "infinity":
                    if max(newRoots) >= hysteresis: # min of the roots is always 0
                        hysteresis = max(newRoots)
                        minBeta = newBeta # Because the epidemic fraction increases with increasing beta
            elif len(newRoots) == 2:
                # testBeta = newBeta - tolerance/2.0
                # testRoots = solveEquilibrium(gamma, testBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, digits=digits) # perturb the solution to see if this is at the end of the branch or not
                # if len(testRoots) == 1:
                #     minBeta = newBeta
                #     # if max(newRoots) >= hysteresis and max(newRoots)>0.1: # min of the roots is always 0
                #     #     hysteresis = max(newRoots)
                # else:
                #     maxBeta = newBeta
                maxBeta = newBeta
            else: # if just the zero solution
                minBeta = newBeta
        return hysteresis
    else:
        return float("nan")


def calculateTheoreticalCriticalAlpha(gamma, minBeta, maxBeta, minAlpha, maxAlpha, degreeHist, meanSimplexDegree=None, isIndependent=False, option="infinity", digits=4, tolerance=0.0001):
    if meanSimplexDegree == None:
        meanSimplexDegree = sum([k*prob for k, prob in degreeHist])

    minAlphaCrit = minAlpha
    maxAlphaCrit = maxAlpha

    minAlphaHysteresis = calculateTheoreticalHysteresis(gamma, minBeta, maxBeta, minAlphaCrit, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, option=option, digits=digits, tolerance=tolerance)

    maxAlphaHysteresis = calculateTheoreticalHysteresis(gamma, minBeta, maxBeta, maxAlphaCrit, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, option=option, digits=digits, tolerance=tolerance)

    if minAlphaHysteresis < tolerance and maxAlphaHysteresis > tolerance:
        # Bisection method
        while maxAlphaCrit - minAlphaCrit > tolerance:
            newAlpha = 0.5*(minAlphaCrit + maxAlphaCrit)
            newHysteresis = calculateTheoreticalHysteresis(gamma, minBeta, maxBeta, newAlpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, option=option, digits=digits, tolerance=0.1*tolerance)
            if newHysteresis == 0:
                minAlphaCrit = newAlpha
            else:
                maxAlphaCrit = newAlpha

        print(0.5*(minAlphaCrit + maxAlphaCrit), flush=True)
        return 0.5*(minAlphaCrit + maxAlphaCrit)
        #return minAlphaCrit
    else:
        print("NaN", flush=True)
        return float("nan")


### Additional functions

# def findHysteresisBoundaries(gamma, betaList, alphaList, minDegree, maxDegree, meanSimplexDegree, tolerance, isIndependent=False, type="uniform", r=4):
#     # These depend on there being an odd number of points
#     n = len(alphaList)
#     m = len(betaList)
#     leftBoundary = []
#     rightBoundary = []
#     alphaHys = []
#     for i in range(n):
#         left = max(betaList)
#         right = 0
#         isHysteresis = False
#         for j in range(m):
#             roots = solveEquilibrium(gamma, betaList[j], alphaList[i], minDegree, maxDegree, meanSimplexDegree, digits=4, isIndependent=isIndependent, type=type, r=r)
#             if isIndependent:
#                 if len(roots) > 2:
#                     isHysteresis = True
#                     if betaList[j] < left:
#                         left = betaList[j]
#
#                     if betaList[j] > right:
#                         right = betaList[j]
#             else:
#                 if len(roots) > 1:
#                     isHysteresis = True
#                     if betaList[j] < left:
#                         left = betaList[j]
#
#                     if betaList[j] > right:
#                         right = betaList[j]
#
#         if isHysteresis:
#             leftBoundary.append(left)
#             rightBoundary.append(right)
#             alphaHys.append(alphaList[i])
#
#     # if type == "uniform":
#     #     meanDegree = 0.5*(minDegree + maxDegree)
#     #
#     # elif type == "power-law":
#     #     meanDegree = avgOfPowerLaw(minDegree, maxDegree, r)
#     #
#     # if not isIndependent:
#     #     meanSimplexDegree = meanDegree
#
#     return leftBoundary, rightBoundary, alphaHys

# def calculateTheoreticalHysteresisOriginal(gamma, betaTheory, alpha, degreeHist, meanSimplexDegree=None, isIndependent=False, option="infinity", digits=4):
#     if meanSimplexDegree == None:
#         meanSimplexDegree = sum([k*prob for k, prob in degreeHist])
#     hysteresis = 0
#     for betaVal in betaTheory:
#         roots = solveEquilibrium(gamma, betaVal, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, digits=digits)
#         if len(roots) == 3:
#             if option == "infinity":
#                 if max(roots)-min(roots) > hysteresis:
#                     hysteresis = max(roots)-min(roots)
#     return hysteresis

# def calculateTheoreticalCriticalAlphaOriginal(gamma, beta, minAlpha, maxAlpha, degreeHist, meanSimplexDegree=None, isIndependent=False, option="infinity", digits=4, tolerance=0.0001):
#     if meanSimplexDegree == None:
#         meanSimplexDegree = sum([k*prob for k, prob in degreeHist])
#
#     minAlphaCrit = minAlpha
#     maxAlphaCrit = maxAlpha
#
#     minAlphaHysteresis = calculateTheoreticalHysteresisOriginal(gamma, beta, minAlphaCrit, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, option=option, digits=digits)
#     maxAlphaHysteresis = calculateTheoreticalHysteresisOriginal(gamma, beta, maxAlphaCrit, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, option=option, digits=digits)
#
#     if minAlphaHysteresis < 0.01 and maxAlphaHysteresis > 0.01:
#         # Bisection method
#         while maxAlphaCrit - minAlphaCrit > tolerance:
#             newAlpha = 0.5*(minAlphaCrit + maxAlphaCrit)
#             newHysteresis = calculateTheoreticalHysteresisOriginal(gamma, beta, newAlpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, option=option, digits=digits)
#             if newHysteresis < 0.01:
#                 minAlphaCrit = newAlpha
#             else:
#                 maxAlphaCrit = newAlpha
#
#         print(0.5*(minAlphaCrit + maxAlphaCrit), flush=True)
#         return 0.5*(minAlphaCrit + maxAlphaCrit)
#     else:
#         print("NaN", flush=True)
#         return float("nan")

#
# def solveEquilibriumOld(gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, digits=4, isIndependent=False, type="power-law", r=4):
#     if isIndependent:
#         return solveIndependentEquilbriumOld(gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, r=r, digits=digits, type=type)
#         #return solveIndependentEquilbriumExpansion(gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, r=r, digits=digits, type=type)
#
#     else:
#         print("Test")
#         return solveDependentEquilbriumOld(gamma, beta, alpha, minDegree, maxDegree, r=r, digits=digits, type=type)
#
#
#
# def solveDependentEquilbriumOld(gamma, beta, alpha, minDegree, maxDegree, r=4, digits=4, type="power-law"):
#     initialGuesses = np.linspace(0, 1, 5)
#     roots = list()
#     for initialGuess in initialGuesses:
#         root, data, ier, msg = fsolve(dependentEquilibriumFunctionOld, initialGuess,  args=(gamma, beta, alpha, minDegree, maxDegree, r, type), full_output=True)
#         if ier == 1:
#             avgInfectionPt = round(calculateAvgInfected(np.asscalar(root), gamma, beta, alpha, minDegree, maxDegree), digits)
#             if avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
#                 roots.append(avgInfectionPt)
#     return roots
#
# def solveIndependentEquilbriumOld(gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, r=4, digits=4, type="power-law"):
#     initialGuesses = [[0, 0], [0.001, 0.001], [0.25, 0.25], [0.5, 0.5],[0.75, 0.75], [1, 1]]
#     roots = list()
#     for initialGuess in initialGuesses:
#         result = root(independentEquilibriumFunctionOld, initialGuess,  args=(gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, r, type))
#         avgInfectionPt = round(np.asscalar(result.x[0]), digits)
#         if result.success and avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
#             roots.append(avgInfectionPt)
#     return roots
#
# def dependentEquilibriumFunctionOld(V, gamma, beta, alpha, minDegree, maxDegree, r, networkDist):
#     if networkDist == "power-law":
#         meanDegree = avgOfPowerLaw(minDegree, maxDegree, r)
#         sum = 0
#         for k in range(minDegree, maxDegree+1):
#             sum = sum + truncatedPowerLaw(k, minDegree, maxDegree, r)*k**2*(beta + alpha*V)/(gamma + beta*k*V + alpha*k*V**2)
#         return 1/meanDegree*sum - 1
#
#     elif networkDist == "uniform":
#         meanDegree = 0.5*(minDegree+maxDegree)
#         frac = 1/(maxDegree-minDegree+1)
#         sum = 0
#         for k in range(minDegree, maxDegree+1):
#             sum = sum + k**2*(beta + alpha*V)/(gamma + beta*k*V + alpha*k*V**2)
#         return frac/meanDegree*sum - 1
#
#     else:
#         print("invalid choice")
#         return
#
# def independentEquilibriumFunctionOld(vars, gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, r, networkDist):
#     # Add this calculation
#     U = vars[0]
#     V = vars[1]
#     if networkDist == "power-law":
#         meanDegree = avgOfPowerLaw(minDegree, maxDegree, r)
#         sumV = 0
#         sumU = 0
#         for k in range(minDegree, maxDegree+1):
#             sumU = sumU + truncatedPowerLaw(k, minDegree, maxDegree, r)*(k*beta*V + alpha*meanSimplexDegree*U**2)/(gamma + beta*k*V + alpha*meanSimplexDegree*U**2)
#             sumV = sumV + truncatedPowerLaw(k, minDegree, maxDegree, r)*k*(k*beta*V + alpha*meanSimplexDegree*U**2)/(gamma + beta*k*V + alpha*meanSimplexDegree*U**2)
#
#         return [sumU-U, 1/meanDegree*sumV - V]
#     elif networkDist == "uniform":
#         meanDegree = 0.5*(minDegree + maxDegree)
#         frac = 1/(maxDegree-minDegree+1)
#         sumV = 0
#         sumU = 0
#         for k in range(minDegree, maxDegree+1):
#             sumU = sumU + (k*beta*V + alpha*meanSimplexDegree*U**2)/(gamma + beta*k*V + alpha*meanSimplexDegree*U**2)
#             sumV = sumV + k*(k*beta*V + alpha*meanSimplexDegree*U**2)/(gamma + beta*k*V + alpha*meanSimplexDegree*U**2)
#
#         return [frac*sumU-U, frac/meanDegree*sumV - V]
#     else:
#         print("invalid choice")
#         return
#
# def solveUniformEquilbrium(gamma, beta, alpha, minDegree, maxDegree, digits=4):
#     initialGuesses = np.linspace(0, 1, 5)
#     roots = list()
#     for initialGuess in initialGuesses:
#         root, data, ier, msg = fsolve(uniformEquilibriumFunction, initialGuess,  args=(gamma, beta, alpha, minDegree, maxDegree), full_output=True)
#         if ier == 1:
#             avgInfectionPt = round(calculateAvgInfected(np.asscalar(root), gamma, beta, alpha, minDegree, maxDegree), digits)
#             if avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
#                 roots.append(avgInfectionPt)
#     return roots
#
# def solvePowerLawEquilbrium(gamma, beta, alpha, minDegree, maxDegree, r, digits=4):
#     initialGuesses = np.linspace(0, 1, 5)
#     roots = list()
#     for initialGuess in initialGuesses:
#         root, data, ier, msg = fsolve(powerLawEquilibriumFunction, initialGuess,  args=(gamma, beta, alpha, minDegree, maxDegree, r), full_output=True)
#         if ier == 1:
#             avgInfectionPt = round(calculateAvgInfected(np.asscalar(root), gamma, beta, alpha, minDegree, maxDegree), digits)
#             if avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
#                 roots.append(avgInfectionPt)
#     return roots
#
#
# def uniformEquilibriumFunction(V, gamma, beta, alpha, minDegree, maxDegree):
#     meanDegree = 0.5*(minDegree+maxDegree)
#     frac = 1/(maxDegree-minDegree+1)
#     sum = 0
#     for k in range(minDegree, maxDegree+1):
#         sum = sum + k**2*(beta + alpha*V)/(gamma + beta*k*V + alpha*k*V**2)
#     return frac/meanDegree*sum - 1
#
# def powerLawEquilibriumFunction(V, gamma, beta, alpha, minDegree, maxDegree, r):
#     meanDegree = avgOfPowerLaw(minDegree, maxDegree, r)
#     sum = 0
#     for k in range(minDegree, maxDegree+1):
#         sum = sum + truncatedPowerLaw(k, minDegree, maxDegree, r)*k**2*(beta + alpha*V)/(gamma + beta*k*V + alpha*k*V**2)
#     return 1/meanDegree*sum - 1
#
# def solveIndependentEquilbriumExpansion(gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, r, digits=4, type="power-law"):
#     initialGuesses = np.linspace(0, 1, 5)
#     roots = list()
#     for initialGuess in initialGuesses:
#         root, data, ier, msg = fsolve(independentExpansionEquilibriumFunction, initialGuess,  args=(gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, r, type), full_output=True)
#         if ier == 1:
#             avgInfectionPt = round(calculateAvgInfected(np.asscalar(root), gamma, beta, alpha, minDegree, maxDegree), digits)
#             if avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
#                 roots.append(avgInfectionPt)
#     return roots
#
#
# def independentExpansionEquilibriumFunction(V, gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, r, networkDist):
#     if networkDist == "power-law":
#         meanDegree = avgOfPowerLaw(minDegree, maxDegree, r)
#         meanSquaredDegree = meanPowerOfPowerLaw(minDegree, maxDegree, r, 2)
#         c = meanDegree**4*meanSimplexDegree/meanSquaredDegree**2
#         sum = 0
#         for k in range(minDegree, maxDegree+1):
#             sum = sum + truncatedPowerLaw(k, minDegree, maxDegree, r)*k*(k*beta + alpha*c*V)/(gamma + beta*k*V + alpha*c*V**2)
#         return 1/meanDegree*sum - 1
#
#     elif networkDist == "uniform":
#         meanDegree = 0.5*(minDegree+maxDegree)
#         meanSquaredDegree = 1./3*(maxDegree**2+minDegree*maxDegree+minDegree**2)
#         frac = 1/(maxDegree-minDegree+1)
#         c = meanDegree**4*meanSimplexDegree/meanSquaredDegree**2
#         sum = 0
#         for k in range(minDegree, maxDegree+1):
#             sum = sum + k**2*(beta + alpha*c*V)/(gamma + beta*k*V + alpha*c*V**2)
#         return frac/meanDegree*sum - 1
#
#     else:
#         print("invalid choice")
#         return
