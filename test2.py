import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs
import random
import simplexUtilities
import simplexContagion
import pickle
from datetime import datetime
import time
import multiprocessing as mp
import simplexTheory
import os

n = 10000
minDegree = 50
maxDegree = n

degreeSequenceRandom = simplexUtilities.generatePowerLawDegreeSequence(n, minDegree, maxDegree, 3, isRandom=True)

#randomA = simplexUtilities.generateConfigModelAdjacency(degreeSequenceRandom)

degreeSequenceNonRandom = simplexUtilities.generatePowerLawDegreeSequence(n, minDegree, maxDegree, 3, isRandom=False)

#nonRandomA = simplexUtilities.generateConfigModelAdjacency(degreeSequenceNonRandom)

print(simplexUtilities.meanPowerOfDegree(degreeSequenceRandom, 1))
print(simplexUtilities.meanPowerOfDegree(degreeSequenceRandom, 2))
print(simplexUtilities.meanPowerOfDegree(degreeSequenceRandom, 3))

print(simplexUtilities.meanPowerOfDegree(degreeSequenceNonRandom, 1))
print(simplexUtilities.meanPowerOfDegree(degreeSequenceNonRandom, 2))
print(simplexUtilities.meanPowerOfDegree(degreeSequenceNonRandom, 3))


percs = np.linspace(0,100,10000)
qn_a = np.percentile(degreeSequenceRandom, percs)
qn_b = np.percentile(degreeSequenceNonRandom, percs)

plt.plot(qn_a,qn_b, ls="", marker="o")
plt.xlabel("Random")
plt.ylabel("Non-random")

x = np.linspace(np.min((qn_a.min(),qn_b.min())), np.max((qn_a.max(),qn_b.max())))
plt.plot(x,x, color="k", ls="--")

plt.show()
#
# plt.figure()
# plt.subplot(121)
# plt.spy(randomA.todense())
#
# plt.subplot(122)
# plt.spy(nonRandomA.todense())
#
# plt.show()
# plt.figure()
# plt.semilogy(degreeSequenceRandom)
# plt.semilogy(degreeSequenceNonRandom)

plt.show()
