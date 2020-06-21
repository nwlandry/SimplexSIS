import sys
sys.path.insert(1, 'C:/Users/nicho/Documents/GitHub/SimplexSIS')

import simplexTheory
import simplexVisualize
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *

filenames = ['Paper Figures/bistability_uniform_uncorrelated','Paper Figures/bistability_r=4_uncorrelated','Paper Figures/bistability_r=3_uncorrelated'] #
#filenames = ['Paper Figures/bistability_uniform_degreecorrelated','Paper Figures/bistability_r=4_degreecorrelated','Paper Figures/bistability_r=3_degreecorrelated']
plt.figure()
colors = ['black','blue','green','maroon', 'grey','red']
i = 0
for file in filenames:
    with open(file, 'rb') as file:
        data = pickle.load(file)
    alphaSim = data[0]
    simulationHysteresis = data[1]
    alphaTheory = data[2]
    theoreticalHysteresis = data[3]
    simulationLabel = data[4]
    theoreticalLabel = data[5]
    plt.plot(alphaSim, simulationHysteresis, 'o-', label=simulationLabel, color=colors[i])
    plt.plot(alphaTheory, theoreticalHysteresis, '-', label=theoreticalLabel, color=colors[i])
    i = i + 1

plt.xlabel(r"$\beta_3$", fontsize=18)
plt.ylabel(r"$B(\beta_3)$", fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12, loc="upper left", frameon=False)
plt.xlim([-0.001,.1])
plt.ylim([0,1])
plt.tight_layout()
plt.show()
