import random
import scipy
from scipy.sparse import csr_matrix
import numpy as np
import multiprocessing as mp

def meanPowerOfDegree(k, power):
    return np.asscalar(np.power(k, power).mean())

def invCDF(u, k0, n, exponent):
    return (k0**(1-exponent) + u*(n**(1-exponent) - k0**(1-exponent)))**(1/(1-exponent))

def generatePowerLawDegreeSequence(numPoints, k0, n, gamma):
    k = []

    for i in range(numPoints):
        u = random.uniform(0, 1)
        k.append(round(invCDF(u, k0, n, gamma)))
    return sorted(k)

def generateUniformDegreeSequence(numPoints, min, max):
    k = []

    for i in range(numPoints):
        u = random.randrange(min, max)
        k.append(round(u))
    return k

def generatePoissonDegreeSequence(numPoints, meanDegree):
    return np.random.poisson(lam=meanDegree, size=numPoints).tolist()


def generateChungLuAdjacency(kIn, kOut, n):
    meanDegree = sum(kIn)/len(kIn)
    if len(kIn)!=n or len(kOut)!=n:
        return

    rowIndex = []
    columnIndex = []
    weights = []
    for i in range(n):
        for j in range(n):
            u = random.uniform(0, 1)
            if u < min(kIn[i]*kOut[j]/(n*meanDegree), 1):
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
