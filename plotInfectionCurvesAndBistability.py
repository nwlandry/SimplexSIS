import pickle
import matplotlib.pyplot as plt
import numpy as np
import simplexUtilities
import simplexVisualize


filename = 'equilibriaData_power-law_r=3_dep_nonrandom_degree'
filename = 'Non-Random Degree/equilibriaData_power-law_r=3_dep_nonrandom_degree'
#filename = 'equilibriaData05172020-013422'
#filename = 'equilibriaData03042020-165346'
#filename = 'Poster/uniform_dep'
with open(filename, 'rb') as file:
    data = pickle.load(file)
gamma = data[0]
beta = data[1]
alpha = data[2]
equilibria = data[3]


plt.figure()
for i in range(len(equilibria)):
    plt.plot(beta, equilibria[i], 'o-', label=r"$\alpha=$" + str(round(alpha[i],3)))
#plt.legend(loc='lower left', fontsize=14)
plt.xlabel(r"$\beta$", fontsize=24)
plt.ylabel("Fraction infected", fontsize=24)
plt.ylim([0, 1])
plt.tight_layout()
plt.show()

bistability = []
for i in range(len(alpha)):
    bistability.append(simplexVisualize.calculateBistability(equilibria[i], beta, option='infinity'))

plt.figure()
plt.plot(alpha, bistability, 'o-', label="Hysteresis")
plt.xlabel(r"$\alpha$", fontsize=24)
plt.ylabel("Hysteresis", fontsize=24)
plt.legend(loc="upper left", fontsize=14)
plt.tight_layout()
plt.show()