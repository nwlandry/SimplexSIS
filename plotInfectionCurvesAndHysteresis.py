import pickle
import matplotlib.pyplot as plt
import numpy as np
import simplexUtilities
import visualizeData


filename = 'equilibriaData_power-law_r=3_indep_final'
filename = 'equilibriaData02252020-130441'
with open(filename, 'rb') as file:
    data = pickle.load(file)
gamma = data[0]
beta = data[1]
alpha = data[2]
equilibria = data[3]

plt.figure()
for i in range(len(equilibria)):
    plt.plot(beta, equilibria[i], 'o-', label=r"$\alpha=$" + str(round(alpha[i],3)))
plt.legend(loc='lower right')
plt.xlabel(r"$\beta$")
plt.ylabel("Fraction infected")
plt.ylim([0,1])
plt.show()

hysteresis = []
for i in range(len(alpha)):
    hysteresis.append(visualizeData.calculateHysteresis(equilibria[i], beta, option='infinity'))

plt.figure()
plt.plot(alpha, hysteresis, 'o-', label="Hysteresis")

plt.xlabel(r"$\alpha$")
plt.ylabel("Hysteresis")
plt.legend(loc="lower right")
plt.show()
