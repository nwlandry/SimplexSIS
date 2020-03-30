import simplexTheory
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *

#filenames = ['hysteresis_uniform_indep_nonrandom_degree','hysteresis_r=4_indep_final','hysteresis_r=3_indep_final'] # Heterogeneity (indep)
#filenames = ['hysteresis_uniform_dep_nonrandom_degree','hysteresis_r=4_dep_final','hysteresis_r=3_dep_nonrandom_degree'] # Heterogeneity (dep)
filenames = ['hysteresis_r=3_indep_final','hysteresis_r=3_dep_nonrandom_degree'] # dep vs. indep
#filenames = ['hysteresis_uniform_indep_nonrandom_degree','hysteresis_r=4_indep_final','hysteresis_r=3_indep_final','hysteresis_uniform_dep_nonrandom_degree','hysteresis_r=4_dep_final','hysteresis_r=3_dep_nonrandom_degree']
plt.figure()
colors = ['black','blue','green','maroon', 'grey','red']
i = 0
for file in filenames:
    with open(file, 'rb') as file:
        data = pickle.load(file)
    alpha = data[0]
    simulationHysteresis = data[1]
    theoreticalHysteresis = data[2]
    simulationLabel = data[3]
    theoreticalLabel = data[4]
    plt.plot(alpha, simulationHysteresis, 'o-', label=simulationLabel, color=colors[i])
    plt.plot(alpha, theoreticalHysteresis, '-', label=theoreticalLabel, color=colors[i])
    i = i + 1

plt.xlabel(r"$\alpha$", fontsize=24)
plt.ylabel("Hysteresis (sup-norm)", fontsize=24)
plt.legend(fontsize=12, loc="upper left")
plt.tight_layout()
plt.show()
