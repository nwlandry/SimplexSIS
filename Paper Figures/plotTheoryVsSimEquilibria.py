import sys
sys.path.insert(1, 'C:/Users/nicho/Documents/GitHub/SimplexSIS')

import pickle
import matplotlib.pyplot as plt
import numpy as np
import simplexUtilities
import simplexVisualize
import simplexTheory

filename = 'equilibriaData06052020-183753'
with open(filename, 'rb') as file:
    data = pickle.load(file)
gamma = data[0]
beta = data[1]
alpha = data[2]
equilibria = data[3]
degreeSequence = data[4]
isDegreeCorrelated = data[5]
degreeDistType = data[6]
exponent = data[7]
meanSimplexDegree = data[8]
print(data[6])
print(data[7])
print(data[5])
index1 = 8
index2 = 16
index3 = 24
print(alpha[index1])
print(alpha[index2])
print(alpha[index3])

degreeHist = simplexTheory.degreeSequenceToHist(degreeSequence)

plt.figure(figsize=(5,10))
plt.subplot(311)
ax1 = plt.gca()
ax1.text(0.54, 0.6, "(a)", fontsize=16)
simplexVisualize.plotTheoreticalAndSimInfectionCurves(ax1, equilibria[index1], gamma, beta, alpha[index1], degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, numTheoryPoints=100, digits=5)

plt.subplot(312)
ax2 = plt.gca()
ax2.text(0.54, 0.6, "(b)", fontsize=16)
simplexVisualize.plotTheoreticalAndSimInfectionCurves(ax2, equilibria[index2], gamma, beta, alpha[index2], degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, numTheoryPoints=100, digits=5)

plt.subplot(313)
ax3 = plt.gca()
ax3.text(0.54, 0.6, "(c)", fontsize=16)
simplexVisualize.plotTheoreticalAndSimInfectionCurves(ax3, equilibria[index3], gamma, beta, alpha[index3], degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, numTheoryPoints=400, digits=5)


ax3.set_xlabel(r'$\beta_2/\beta_2^c$', fontsize=16)
ax1.set_ylabel(r'$U$', fontsize=16)
ax2.set_ylabel(r'$U$', fontsize=16)
ax3.set_ylabel(r'$U$', fontsize=16)
plt.tight_layout()
plt.savefig('equilibria_fig.svg',bbox_inches='tight',dpi = 600)
plt.show()
