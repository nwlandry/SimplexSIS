import pickle
import matplotlib.pyplot as plt
import numpy as np
import simplexUtilities
import simplexContagion


filename = 'equilibriaData12262019-000110'
with open(filename, 'rb') as file:
    data = pickle.load(file)
alpha = data[0]
beta = data[1]
gamma = data[2]
kAvg = data[3]
kAvgSimplex = data[4]
equilibria = data[5]
alphaCrit = data[6]
betaCrit = data[7]

plt.figure()
for i in range(len(equilibria)):
    plt.plot(beta, equilibria[i], 'o-', label=r"$\lambda_{\alpha}=$" + str(round(alpha[i],3)))
plt.legend(loc='lower right')
plt.xlabel(r"$\beta$")
plt.ylabel("Fraction infected")
plt.ylim([0,1])
plt.show()

hysteresis = []
for i in range(len(alpha)):
    hysteresis.append(simplexContagion.calculateHysteresis(equilibria[i], beta, option='infinity'))

plt.figure()
plt.plot(alpha, hysteresis, 'o-', label="Hysteresis")
plt.scatter(alphaCrit, 0, s=50, label=r"$\alpha_{crit}$",  facecolors='none', edgecolors='red')
plt.xlabel(r"$\alpha$")
plt.ylabel("Hysteresis")
plt.legend(loc="lower right")
plt.show()
