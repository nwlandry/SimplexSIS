from scipy.optimize import fsolve, root
import matplotlib.pyplot as plt
import numpy as np

def getPhase(gamma, beta, alpha, degreeHist, meanSimplexDegree=None, isDegreeCorrelated=True, majorityVote=True, healing=False, digits=4):
    return len(solveEquilibrium(gamma, beta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, majorityVote=majorityVote, healing=healing, digits=digits))

def solveEquilibrium(gamma, beta, alpha, degreeHist, meanSimplexDegree=None, isDegreeCorrelated=True, majorityVote=True, healing=False, digits=4):
    if meanSimplexDegree is None:
        meanSimplexDegree = generateMeanDegreeFromHist(degreeHist)

    if majorityVote:
        if isDegreeCorrelated:
            return solveDegreeCorrelatedEquilbriumMajorityVote(gamma, beta, alpha, degreeHist, healing=healing, digits=digits)
        else:
            return solveUncorrelatedEquilbriumMajorityVote(gamma, beta, alpha, degreeHist, meanSimplexDegree, healing=healing, digits=digits)

    else:
        if isDegreeCorrelated:
            return solveDegreeCorrelatedEquilbriumIndividual(gamma, beta, alpha, degreeHist, healing=healing, digits=digits)
        else:
            return solveUncorrelatedEquilbriumIndividual(gamma, beta, alpha, degreeHist, meanSimplexDegree, healing=healing, digits=digits)


def solveDegreeCorrelatedEquilbriumMajorityVote(gamma, beta, alpha, degreeHist, healing=False, digits=4):
    initialGuesses = [val**2 for val in np.linspace(0, 0.8, 10)]
    roots = list()
    for initialGuess in initialGuesses:
        root, data, ier, msg = fsolve(degreeCorrelatedEquilibriumFunctionMajorityVote, initialGuess,  args=(gamma, beta, alpha, degreeHist, healing), full_output=True)
        if ier == 1:
            avgInfectionPt = round(calculateMeanInfectedDegreeCorrelatedMajorityVote(np.asscalar(root), gamma, beta, alpha, degreeHist, healing), digits)
            if avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
                roots.append(avgInfectionPt)
    return roots

def solveUncorrelatedEquilbriumMajorityVote(gamma, beta, alpha, degreeHist, meanSimplexDegree, healing=False, digits=4):
    initialGuesses = [[val**2,val**2] for val in np.linspace(0, 0.8, 10)]
    roots = list()
    for initialGuess in initialGuesses:
        result = root(uncorrelatedEquilibriumFunctionMajorityVote, initialGuess,  args=(gamma, beta, alpha, degreeHist, meanSimplexDegree, healing))
        avgInfectionPt = round(np.asscalar(result.x[0]), digits)
        if result.success and avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
            roots.append(avgInfectionPt)
    return roots

def solveDegreeCorrelatedEquilbriumIndividual(gamma, beta, alpha, degreeHist, healing=False, digits=4):
    initialGuesses = [val**2 for val in np.linspace(0, 0.8, 10)]
    roots = list()
    for initialGuess in initialGuesses:
        root, data, ier, msg = fsolve(degreeCorrelatedEquilibriumFunctionIndividual, initialGuess,  args=(gamma, beta, alpha, degreeHist, healing), full_output=True)
        if ier == 1:
            avgInfectionPt = round(calculateMeanInfectedDegreeCorrelatedIndividual(np.asscalar(root), gamma, beta, alpha, degreeHist, healing), digits)
            if avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
                roots.append(avgInfectionPt)
    return roots

def solveUncorrelatedEquilbriumIndividual(gamma, beta, alpha, degreeHist, meanSimplexDegree, healing=False, digits=4):
    initialGuesses = [[val**2,val**2] for val in np.linspace(0, 0.8, 10)]
    roots = list()
    for initialGuess in initialGuesses:
        result = root(uncorrelatedEquilibriumFunctionIndividual, initialGuess,  args=(gamma, beta, alpha, degreeHist, meanSimplexDegree, healing))
        avgInfectionPt = round(np.asscalar(result.x[0]), digits)
        if result.success and avgInfectionPt not in set(roots) and avgInfectionPt <= 1 and avgInfectionPt >= 0:
            roots.append(avgInfectionPt)
    return roots

def degreeCorrelatedEquilibriumFunctionMajorityVote(V, gamma, beta, alpha, degreeHist, healing=False):
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

def uncorrelatedEquilibriumFunctionMajorityVote(vars, gamma, beta, alpha, degreeHist, meanSimplexDegree, healing=False):
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

def degreeCorrelatedEquilibriumFunctionIndividual(V, gamma, beta, alpha, degreeHist, healing=False):
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

def uncorrelatedEquilibriumFunctionIndividual(vars, gamma, beta, alpha, degreeHist, meanSimplexDegree, healing=False):
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

def generateTheoreticalDegreeHist(minDegree, maxDegree, networkDist, exponent=4):
    degreeHist = list()
    for degree in range(minDegree, maxDegree+1):
        if networkDist == "power-law":
            prob = truncatedPowerLaw(degree, minDegree, maxDegree, exponent)
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

def calculateMeanInfectedDegreeCorrelatedMajorityVote(V, gamma, beta, alpha, degreeHist, healing=False):
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

def calculateMeanInfectedDegreeCorrelatedIndividual(V, gamma, beta, alpha, degreeHist, healing=False):
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

def calculateMeanInfectedUncorrelatedFromV(V, gamma, beta, alpha, degreeHist, meanSimplexDegree):
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

def truncatedPowerLaw(k, minDegree, maxDegree, exponent):
    return (exponent-1)/(minDegree**(1-exponent)-maxDegree**(1-exponent))*k**(-exponent)

def avgOfPowerLaw(minDegree, maxDegree, exponent):
    return (minDegree**(2-exponent)-maxDegree**(2-exponent))*(exponent-1)/((minDegree**(1-exponent)-maxDegree**(1-exponent))*(exponent-2))

def avgOfPowerLawEqn(minDegree, maxDegree, exponent, meanDeg):
    return (minDegree**(2-exponent)-maxDegree**(2-exponent))*(exponent-1)/((minDegree**(1-exponent)-maxDegree**(1-exponent))*(exponent-2)) - meanDeg

def meanPowerOfPowerLaw(minDegree, maxDegree, exponent, power):
    return (minDegree**(power+1-exponent)-maxDegree**(power+1-exponent))*(pwoer-1)/((minDegree**(1-exponent)-maxDegree**(1-exponent))*(exponent-power-1))

def calculateTheoreticalBistabilityVisually(gamma, minBeta, maxBeta, alpha, degreeHist, meanSimplexDegree=None, isDegreeCorrelated=True, digits=4, tolerance=0.0001, stopAtBistability=False):
    if meanSimplexDegree == None:
        meanSimplexDegree = sum([k*prob for k, prob in degreeHist])

    plt.figure()
    minRoots = solveEquilibrium(gamma, minBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits)

    maxRoots = solveEquilibrium(gamma, maxBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits)

    plt.scatter(minBeta, len(minRoots))
    plt.scatter(maxBeta, len(maxRoots))
    bistabilityIndex = 0
    if (max(minRoots) < tolerance and max(maxRoots) > tolerance) or len(minRoots) == 3:
        iter = 0
        while maxBeta - minBeta > tolerance:
            newBeta = 0.5*(minBeta + maxBeta)
            newRoots = solveEquilibrium(gamma, newBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits)
            plt.scatter(newBeta,len(newRoots))
            if len(newRoots) == 3:
                if option == "infinity":
                    if max(newRoots) >= bistabilityIndex: # min of the roots is always 0
                        bistabilityIndex = max(newRoots)
                        if stopAtBistability:
                            break
                        minBeta = newBeta # Because the epidemic fraction increases with increasing beta
            elif len(newRoots) == 2:
                testBeta = newBeta - tolerance/2.0
                testRoots = solveEquilibrium(gamma, testBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits) # perturb the solution to see if this is at the end of the branch or not
                if len(testRoots) == 1:
                    minBeta = newBeta
                else:
                    maxBeta = newBeta
            else: # if just the zero solution
                minBeta = newBeta
            iter = iter + 1
        print(iter)
        plt.show()
        return bistabilityIndex
    else:
        return float("nan")

def calculateTheoreticalBistability(gamma, betaCrit, alpha, degreeHist, meanSimplexDegree=None, isDegreeCorrelated=True, digits=4, tolerance=0.0001, stopAtBistability=False, option="fast"):
    if meanSimplexDegree == None:
        meanSimplexDegree = sum([k*prob for k, prob in degreeHist])

    minBeta = 0.5*betaCrit
    maxBeta = 1.5*betaCrit

    if option =="fast":
        roots = solveEquilibrium(gamma, betaCrit-tolerance, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits)
        if len(roots) == 3:
            return max(roots)
        else:
            return 0

    elif option == "bisection":
        minRoots = solveEquilibrium(gamma, minBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits)

        maxRoots = solveEquilibrium(gamma, maxBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits)
        bistabilityIndex = 0
        if (max(minRoots) < tolerance and max(maxRoots) > tolerance) or len(minRoots) == 3:
            while maxBeta - minBeta > tolerance:
                newBeta = 0.5*(minBeta + maxBeta)
                newRoots = solveEquilibrium(gamma, newBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits)

                if len(newRoots) == 3:
                    if max(newRoots) >= bistabilityIndex: # min of the roots is always 0
                        bistabilityIndex = max(newRoots)
                        if stopAtBistability: # Stop if there is bistability to save time
                            break
                        minBeta = newBeta # Because the epidemic fraction increases with increasing beta
                elif len(newRoots) == 2:
                    # taking care of the singularity case
                    testBeta = newBeta + tolerance
                    testRoots = solveEquilibrium(gamma, testBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits)
                    if len(testRoots) == 3:
                        minBeta = newBeta
                        if stopAtBistability:
                            bistabilityIndex = max(testRoots)
                            break
                    else:
                        maxBeta = newBeta
                else: # if just the zero solution
                    minBeta = newBeta
            return bistabilityIndex
        else:
            return float("nan")

    elif option == "trisection":
        minRoots = solveEquilibrium(gamma, minBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits)

        maxRoots = solveEquilibrium(gamma, maxBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits)
        bistabilityIndex = 0

        if len(minRoots) == 3:
            bistabilityIndex = max(minRoots)
        elif (max(minRoots) < tolerance and max(maxRoots) > tolerance):
            bistabilityIndex = 0
        else:
            return float("nan")

        while maxBeta - minBeta > tolerance:
            newMinBeta = 2/3*minBeta + 1/3*maxBeta
            newMaxBeta = 1/3*minBeta + 2/3*maxBeta
            newMinRoots = solveEquilibrium(gamma, newMinBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits)
            newMaxRoots = solveEquilibrium(gamma, newMaxBeta, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits)


            if len(newMinRoots) == 1 and len(newMaxRoots) == 1:
                minBeta = newMaxBeta
            elif len(newMinRoots) == 1 and len(newMaxRoots) == 2:
                minBeta = newMinBeta
            elif len(newMinRoots) == 1 and len(newMaxRoots) == 3:
                if stopAtBistability:
                    bistabilityIndex = max(newMaxRoots)
                    break
                minBeta = newMaxBeta
            elif len(newMinRoots) == 2 and len(newMaxRoots) == 2:
                maxBeta = newMaxBeta
            elif len(newMinRoots) == 2 and len(newMaxRoots) == 3:
                bistabilityIndex = max(newMaxRoots)
                if stopAtBistability:
                    break
                minBeta = newMaxBeta
            elif len(newMinRoots) == 3 and len(newMaxRoots) == 2:
                bistabilityIndex = max(newMinRoots)
                if stopAtBistability:
                    break
                minBeta = newMinBeta
                maxBeta = newMaxBeta
            elif len(newMinRoots) == 3 and len(newMaxRoots) == 3:
                bistabilityIndex = max(newMaxRoots)
                if stopAtBistability:
                    break
                minBeta = newMaxBeta

        return bistabilityIndex

    else:
        print("Invalid choice")

def calculateTheoreticalCriticalAlpha(gamma, betaCrit, minAlpha, maxAlpha, degreeHist, meanSimplexDegree=None, isDegreeCorrelated=True, digits=4, tolerance=0.0001, option="fast"):
    if meanSimplexDegree == None:
        meanSimplexDegree = sum([k*prob for k, prob in degreeHist])

    minAlphaCrit = minAlpha
    maxAlphaCrit = maxAlpha

    bistabilityOfMinAlpha = calculateTheoreticalBistability(gamma, betaCrit, minAlphaCrit, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits, tolerance=tolerance, stopAtBistability=False, option=option)

    bistabilityOfMaxAlpha = calculateTheoreticalBistability(gamma, betaCrit, maxAlphaCrit, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits, tolerance=tolerance, stopAtBistability=False, option=option)

    if bistabilityOfMinAlpha < tolerance and bistabilityOfMaxAlpha > tolerance:
        # Bisection method
        while maxAlphaCrit - minAlphaCrit > tolerance:
            newAlpha = 0.5*(minAlphaCrit + maxAlphaCrit)
            bistabilityOfNewAlpha = calculateTheoreticalBistability(gamma, betaCrit, newAlpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits, tolerance=0.1*tolerance, stopAtBistability=False, option=option)
            if bistabilityOfNewAlpha == 0:
                minAlphaCrit = newAlpha
            else:
                maxAlphaCrit = newAlpha

        print(minAlphaCrit, flush=True) # this is the largest value without bistability
        return minAlphaCrit
    else:
        print("NaN", flush=True)
        return float("nan")
