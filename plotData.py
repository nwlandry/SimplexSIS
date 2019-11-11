import pickle
import matplotlib.pyplot as plt
import numpy as np
import simplexUtilities
import simplexContagion


filename = 'equilibriaData11112019-053507'
with open(filename, 'rb') as file:
    data = pickle.load(file)
alpha = data[0]
beta = data[1]
equilibria = data[2]
plt.figure()
for i in range(len(equilibria)):
    plt.plot(beta, equilibria[i], 'o-', label=r"$\alpha=$" + str(round(alpha[i],3)))
plt.legend()
plt.xlabel(r"$\beta$")
plt.ylabel("Fraction infected")
plt.ylim([0,1])
plt.show()

hysteresis = []
for i in range(len(alpha)):
    hysteresis.append(simplexContagion.calculateHysteresis(equilibria[i], beta, option='area'))

plt.figure()
plt.plot(alpha, hysteresis, 'o-')
plt.xlabel(r"$\alpha$")
plt.ylabel("Hysteresis")
plt.show()
