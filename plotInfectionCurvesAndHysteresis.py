import pickle
import matplotlib.pyplot as plt
import numpy as np
import simplexUtilities
import visualizeData


filename = 'equilibriaData_power-law_r=3_indep_final'
#filename = 'equilibriaData02222020-191255#'
#filename = 'equilibriaData02292020-015156'
#filename = 'Poster/uniform_dep'
with open(filename, 'rb') as file:
    data = pickle.load(file)
gamma = data[0]
beta = data[1]
alpha = data[2]
equilibria = data[3]
# print(data)
# alpha = data[0]
# beta = data[1]
# equilibria = data[5]

plt.figure()
for i in range(len(equilibria)):
    plt.plot(beta, equilibria[i], 'o-', label=r"$\alpha=$" + str(round(alpha[i],3)))
plt.legend(loc='lower left', fontsize=14)
plt.xlabel(r"$\beta$", fontsize=18)
plt.ylabel("Fraction infected", fontsize=18)
plt.ylim([0, 1])
plt.show()

hysteresis = []
for i in range(len(alpha)):
    hysteresis.append(visualizeData.calculateHysteresis(equilibria[i], beta, option='infinity'))

plt.figure()
plt.plot(alpha, hysteresis, 'o-', label="Hysteresis")
plt.xlabel(r"$\alpha$", fontsize=18)
plt.ylabel("Hysteresis", fontsize=18)
plt.legend(loc="upper left", fontsize=14)
plt.show()
