import pickle
import matplotlib.pyplot as plt
import numpy as np
import simplexUtilities
import simplexVisualize


#filename = 'equilibriaData_power-law_r=3_dep_nonrandom_degree'
#filename = 'Non-Random Degree/equilibriaData_power-law_r=3_dep_nonrandom_degree'
#filename = 'equilibriaData06062020-002032'# Uniform DC
#filename = 'equilibriaData06052020-221523' # Uniform UC
#filename = 'equilibriaData06052020-204903' # r=4 DC
#filename = 'equilibriaData06062020-003754' # r=3 DC
filename = 'equilibriaData06112020-050847' # r=3 UC
filename = 'equilibriaData06112020-074534' # r=4 UC
filename = 'equilibriaData06142020-212240' # r=4 UC
with open(filename, 'rb') as file:
    data = pickle.load(file)
gamma = data[0]
beta = data[1]
alpha = data[2]
equilibria = data[3]
degreeSequence = data[4]
print(data[6])
print(data[7])
print(data[5])

plt.figure()
for i in range(len(equilibria)):
    plt.plot(beta, equilibria[i], 'o-', label=r"$\alpha=$" + str(round(alpha[i],3)))
#plt.legend(loc='lower left', fontsize=14)
plt.xlabel(r"$\beta$", fontsize=24)
plt.ylabel("Fraction infected", fontsize=24)
#plt.ylim([0, 1])
plt.tight_layout()
plt.show()

bistability = []
for i in range(len(alpha)):
    bistability.append(simplexVisualize.calculateBistability(equilibria[i], beta, option='infinity'))

plt.figure()
plt.plot(alpha, bistability, 'o-')
plt.xlabel(r"$\beta_3$", fontsize=24)
plt.ylabel(r"$B(\beta_3)$", fontsize=24)
plt.tight_layout()
plt.show()
