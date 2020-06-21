import sys
sys.path.insert(1, 'C:/Users/nicho/Documents/GitHub/SimplexSIS')

import pickle
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mtick

filename = "Paper Figures/phasePlotIndividualContagionHealingFinal"
with open(filename, 'rb') as file:
    data = pickle.load(file)

xMin = data[0]
xMax = data[1]
yMin = data[2]
yMax = data[3]
phaseGrid = data[4]
print(xMin)
print(yMin)
fig = plt.figure()
ax = fig.add_subplot(111)
c = plt.imshow(np.flipud(phaseGrid), interpolation="none", cmap="Reds", extent=[xMin, xMax, yMin, yMax], aspect="auto", vmin=1, vmax=3)
ax.set_xlabel(r"$\beta_2$", fontsize=18)
ax.set_ylabel(r"$\beta_3$", fontsize=18)
plt.xticks([xMin, xMax], fontsize=14)
plt.yticks([yMin, yMax], fontsize=14)
ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))
plt.tight_layout()
plt.show()
