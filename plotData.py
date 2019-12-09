import pickle
import matplotlib.pyplot as plt
import numpy as np
import simplexUtilities
import simplexContagion


filename = 'equilibriaData12042019-233651'
with open(filename, 'rb') as file:
    data = pickle.load(file)
alpha = data[0]
beta = data[1]
gamma = data[2]
kAvg = data[3]
kAvgSimplex = data[4]
equilibria = data[5]

lambdaNetwork = [b*kAvg/gamma for b in beta]
lambdaSimplex = [a*kAvgSimplex/gamma for a in alpha]
plt.figure()
for i in range(len(equilibria)):
    plt.plot(lambdaNetwork, equilibria[i], 'o-', label=r"$\lambda_{\alpha}=$" + str(round(lambdaSimplex[i],3)))
plt.legend(loc='lower left')
plt.xlabel(r"$\lambda$")
plt.ylabel("Fraction infected")
plt.ylim([0,1])
plt.show()

hysteresis = []
for i in range(len(lambdaSimplex)):
    hysteresis.append(simplexContagion.calculateHysteresis(equilibria[i], lambdaNetwork, option='infinity'))

plt.figure()
plt.plot(lambdaSimplex, hysteresis, 'o-')
plt.xlabel(r"$\lambda_{\alpha}$")
plt.ylabel("Hysteresis")
plt.show()
