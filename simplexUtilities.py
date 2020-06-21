import random
import scipy
from simplexTheory import *
from scipy.sparse import csr_matrix
import numpy as np
import multiprocessing as mp

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

        def f(beta3, k1, k2, k3, k4, k5, gamma):
            return a3(beta3, k1, k2, k3, k4, gamma)**2 - 4*a2(beta3, k1, k2, k3, gamma)*a4(beta3, k1, k2, k3, k4, k5, gamma)

        b1 = fsolve(f, 0.1, args=(k1, k2, k3, k4, k5, gamma))

        # saddle-node bifurcation
        if a4(b1, k1, k2, k3, k4, k5, gamma) < 0 and -a3(b1, k1, k2, k3, k4, gamma)/(2*a4(b1, k1, k2, k3, k4, k5, gamma)) > 0 and -a3(b1, k1, k2, k3, k4, gamma)/(2*a4(b1, k1, k2, k3, k4, k5, gamma)) < 1:
            return b1
        else:
            return k3/k1**4*gamma

def meanPowerOfDegree(degreeSequence, power):
    return np.asscalar(np.power(degreeSequence, power).mean())


def invCDFPowerLaw(u, minDegree, maxDegree, exponent):
    return (minDegree**(1-exponent) + u*(maxDegree**(1-exponent) - minDegree**(1-exponent)))**(1/(1-exponent))


def generatePowerLawDegreeSequence(numPoints, minDegree, maxDegree, exponent, isRandom=True):
    degreeSequence = list()
    if isRandom:
        for i in range(numPoints):
            u = random.uniform(0, 1)
            degreeSequence.append(round(invCDFPowerLaw(u, minDegree, maxDegree, exponent)))
        return sorted(degreeSequence)
    else:
        for i in range(numPoints):
            degreeSequence.append(round(invCDFPowerLaw(i/(numPoints-1), minDegree, maxDegree, exponent)))
        return degreeSequence


def generateUniformDegreeSequence(numPoints, minDegree , maxDegree, isRandom=True):
    degreeSequence = list()
    if isRandom:
        for i in range(numPoints):
            u = random.randrange(round(minDegree), round(maxDegree))
            degreeSequence.append(round(u))
        return sorted(degreeSequence)
    else:
        for i in range(numPoints):
            degreeSequence.append(round(minDegree + i/(numPoints-1)*(maxDegree-minDegree)))
        return degreeSequence

def generatePoissonDegreeSequence(numPoints, meanDegree):
    return np.random.poisson(lam=meanDegree, size=numPoints).tolist()


def generateChungLuAdjacency(inDegree, outDegree, n):
    meanDegree = sum(inDegree)/len(inDegree)
    if len(inDegree)!=n or len(outDegree)!=n:
        return

    rowIndex = []
    columnIndex = []
    weights = []
    for i in range(n):
        for j in range(n):
            u = random.uniform(0, 1)
            if u < min(inDegree[i]*outDegree[j]/(n*meanDegree), 1):
                rowIndex.append(i)
                columnIndex.append(j)
                weights.append(1)
    return csr_matrix((np.asarray(weights), (np.asarray(rowIndex), np.asarray(columnIndex))), shape=(n, n))


def generateConfigModelAdjacency(degreeSequence):
    k = degreeSequence.copy()
    n = len(degreeSequence)
    if (sum(k) % 2) != 0:
        i = random.randrange(len(k))
        k[i] = k[i] + 1
    stubs = []
    for index in range(len(k)):
        stubs.extend([index]*int(k[index]))
    rowIndex = []
    columnIndex = []
    weights = []

    while len(stubs) != 0:
        u = random.sample(range(len(stubs)), 2)
        i = stubs[u[0]]
        j = stubs[u[1]]
        rowIndex.extend([i, j])
        columnIndex.extend([j, i])
        # By default csr_matrix adds the duplicates together to account for multiedges
        weights.extend([1, 1])
        for index in sorted(u, reverse=True):
            del stubs[index]
    return csr_matrix((np.asarray(weights), (np.asarray(rowIndex), np.asarray(columnIndex))), shape=(n, n))


def generateConfigModelSimplexList(degreeSequence, simplexSize):
    k = degreeSequence.copy()
    # Making sure we have the right number of stubs
    if (sum(k) % simplexSize) != 0:
        remainder = sum(k) % simplexSize
        for i in range(int(round(simplexSize - remainder))):
            j = random.randrange(len(k))
            k[j] = k[j] + 1

    stubs = []
    simplexList = []
    simplexIndices = [[] for i in range(len(k))]
    # Creating the list to index through
    for index in range(len(k)):
        stubs.extend([index]*int(k[index]))

    simplexIndex = 0
    while len(stubs) != 0:
        u = random.sample(range(len(stubs)), simplexSize)
        simplex = []
        for index in u:
            simplex.append(stubs[index])
            simplexIndices[stubs[index]].append(simplexIndex)

        simplexIndex = simplexIndex + 1
        simplexList.append(simplex)

        for index in sorted(u, reverse=True):
            del stubs[index]

    return simplexList, simplexIndices


def generateUniformSimplexList(n, meanSimplexDegree, simplexSize):
    numSimplices = int(meanSimplexDegree*n/simplexSize)
    simplexIndices = [[] for i in range(n)]
    simplexList = []
    for i in range(numSimplices):
        u = random.choices(range(n), k=simplexSize)
        for index in u:
            simplexIndices[index].append(i)
        simplexList.append(u)
        # I think the chances are very low to have a duplicate.
    return simplexList, simplexIndices
