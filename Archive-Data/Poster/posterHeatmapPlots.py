import simplexTheory
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np


filename = 'Poster/uniform_indep'
#filename = 'Poster/power-law_r=4_indep'
#filename = 'Poster/power-law_r=4point5_indep'
with open(filename, 'rb') as file:
    data1 = pickle.load(file)
alpha1 = data1[0]
beta1 = data1[1]
gamma1 = data1[2]
kAvg1 = data1[3]
kAvgSimplex1 = data1[4]
equilibria1 = data1[5]
alphaCrit1 = data1[6]
betaCrit1 = data1[7]

filename = 'Poster/uniform_dep'
#filename = 'Poster/power-law_r=4_dep'
#filename = 'Poster/power-law_r=4point5_dep'
with open(filename, 'rb') as file:
    data2 = pickle.load(file)
alpha2 = data2[0]
beta2 = data2[1]
gamma2 = data2[2]
kAvg2 = data2[3]
kAvgSimplex2 = data2[4]
equilibria2 = data2[5]
alphaCrit2 = data2[6]
betaCrit2 = data2[7]

minK1 = 14
maxK1 = 100
isIndependent1 = True
type1 = "power-law"
meanSimplexDegree1 = 20
r1 = 4.5

minK2 = 14
maxK2 = 100
isIndependent2 = False
type2 = "power-law"
meanSimplexDegree2 = 20
r2 = 4.5
tolerance = 1e-5

# beta2 = np.linspace(0,0.3,100)
# alpha2 = np.linspace(0,0.3,100)
# plt.figure()
# for i in range(len(beta2)):
#     roots = simplexTheory.solveEquilibrium(gamma2, beta2[i], 0, minK2, maxK2, meanSimplexDegree2, digits=4, isIndependent=isIndependent2, type=type2, r=r2)
#     for root in roots:
#         plt.plot(beta2[i], root, 'k.')
#
# plt.show()

diff1 = simplexContagion.findDiffMatrix(equilibria1)
diff2 = simplexContagion.findDiffMatrix(equilibria2)
# betaTheory = np.linspace(0, 0.1, 40)
# alphaTheory = np.linspace(0.1, 0.264, 20)
# leftBoundary1, rightBoundary1, alphaHys1 = simplexTheory.findHysteresisBoundaries(gamma1, betaTheory, alphaTheory, minK1, maxK1, meanSimplexDegree1, tolerance, isIndependent=isIndependent1, type=type1, r=r1)
# leftBoundary2, rightBoundary2, alphaHys2 = simplexTheory.findHysteresisBoundaries(gamma2, betaTheory, alphaTheory, minK2, maxK2, meanSimplexDegree2, tolerance, isIndependent=isIndependent2, type=type2, r=r2)

#visualizeData.plotTheoryAndSim(alpha2, beta2, diff2, leftBoundary2, rightBoundary2, alphaHys2)
#visualizeData.plotTheoryAndSimCompare(alpha1, beta1, diff1, leftBoundary1, rightBoundary1, alphaHys1, alpha2, beta2, diff2, leftBoundary2, rightBoundary2, alphaHys2)

visualizeData.plotCritPointsAndSimCompare(alpha1, beta1, diff1, alphaCrit1, betaCrit1, alpha2, beta2, diff2, alphaCrit2, betaCrit2)
