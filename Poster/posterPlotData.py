import pickle
import matplotlib.pyplot as plt
import numpy as np
import simplexUtilities
import simplexContagion


filenameDep1 = 'Poster/power-law_r=4_dep'
with open(filenameDep1, 'rb') as file:
    dataDep1 = pickle.load(file)
alphaDep1 = dataDep1[0]
betaDep1 = dataDep1[1]
gammaDep1 = dataDep1[2]
kAvgDep1 = dataDep1[3]
kAvgSimplexDep1 = dataDep1[4]
equilibriaDep1 = dataDep1[5]
alphaCritDep1 = dataDep1[6]
betaCritDep1 = dataDep1[7]

filenameIndep1 = 'Poster/power-law_r=4_indep'
with open(filenameIndep1, 'rb') as file:
    dataIndep1 = pickle.load(file)
alphaIndep1 = dataIndep1[0]
betaIndep1 = dataIndep1[1]
gammaIndep1 = dataIndep1[2]
kAvgIndep1 = dataIndep1[3]
kAvgSimplexIndep1 = dataIndep1[4]
equilibriaIndep1 = dataIndep1[5]
alphaCritIndep1 = dataIndep1[6]
betaCritIndep1 = dataIndep1[7]

filenameDep2 = 'Poster/uniform_dep'
filenameDep2 = 'Poster/uniform_dep_refined'
with open(filenameDep2, 'rb') as file:
    dataDep2 = pickle.load(file)
alphaDep2 = dataDep2[0]
betaDep2 = dataDep2[1]
gammaDep2 = dataDep2[2]
kAvgDep2 = dataDep2[3]
kAvgSimplexDep2 = dataDep2[4]
equilibriaDep2 = dataDep2[5]
alphaCritDep2 = dataDep2[6]
betaCritDep2 = dataDep2[7]

filenameIndep2 = 'Poster/uniform_indep'
filenameIndep2 = 'Poster/uniform_indep_refined'
with open(filenameIndep2, 'rb') as file:
    dataIndep2 = pickle.load(file)
alphaIndep2 = dataIndep2[0]
betaIndep2 = dataIndep2[1]
gammaIndep2 = dataIndep2[2]
kAvgIndep2 = dataIndep2[3]
kAvgSimplexIndep2 = dataIndep2[4]
equilibriaIndep2 = dataIndep2[5]
alphaCritIndep2 = dataIndep2[6]
betaCritIndep2 = dataIndep2[7]

plt.figure()
#indices = range(len(equilibria))
indices = [0, 11, 23]
c = ["green", "blue", "red"]
for i in range(len(indices)):
    plt.plot(betaDep2, equilibriaDep2[indices[i]], '-', linewidth=2, label=r"$\alpha=$" + str(round(alphaDep2[indices[i]],3)), color=c[i])
plt.legend(loc='lower left')
plt.xlabel(r"$\beta$")
plt.ylabel("Fraction infected")
plt.ylim([0,0.55])
plt.show()

hysteresisDep1 = []
hysteresisIndep1 = []
hysteresisDep2 = []
hysteresisIndep2 = []
for i in range(len(alphaDep1)):
    hysteresisDep1.append(simplexContagion.calculateHysteresis(equilibriaDep1[i], betaDep1, option='infinity'))

for i in range(len(alphaIndep1)):
    hysteresisIndep1.append(simplexContagion.calculateHysteresis(equilibriaIndep1[i], betaIndep1, option='infinity'))

for i in range(len(alphaDep2)):
    hysteresisDep2.append(simplexContagion.calculateHysteresis(equilibriaDep2[i], betaDep2, option='infinity'))

for i in range(len(alphaIndep2)):
    hysteresisIndep2.append(simplexContagion.calculateHysteresis(equilibriaIndep2[i], betaIndep2, option='infinity'))

plt.figure()
plt.plot(alphaDep1, hysteresisDep1, '-', linewidth=2, label="r=4 (dep.)", color="red")
plt.plot(alphaIndep1, hysteresisIndep1, '-', linewidth=2, label="r=4 (indep.)", color="green")
plt.plot(alphaDep2, hysteresisDep2, '-', linewidth=2, label="uniform (dep.)", color="blue")
plt.plot(alphaIndep2, hysteresisIndep2, '-', linewidth=2, label="uniform (indep.)", color="orange")

plt.scatter(alphaCritDep1, 0.02, linewidth=2, s=100, marker='o', label=r"$\alpha_{crit}$ r=4 (dep.)", facecolors='none', edgecolors='red')
plt.scatter(alphaCritIndep1, 0.02, s=100, linewidth=2, marker='o', label=r"$\alpha_{crit}$ r=4 (indep.)", facecolors='none', edgecolors='green')
plt.scatter(alphaCritDep2, 0.02, s=100, linewidth=2, marker='o', label=r"$\alpha_{crit}$ uniform (dep.)", facecolors='none', edgecolors='blue')
plt.scatter(alphaCritIndep2, 0.02, s=100, linewidth=2, marker='o', label=r"$\alpha_{crit}$ uniform (indep.)", facecolors='none', edgecolors='orange')
plt.xlim([0, 0.3])
plt.ylim([0, 0.55])
plt.xlabel(r"$\alpha$")
plt.ylabel("Hysteresis")
plt.legend(loc="upper left")
plt.show()
