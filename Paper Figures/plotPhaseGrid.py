import sys
sys.path.insert(1, 'C:/Users/nicho/Documents/GitHub/SimplexSIS')

import pickle
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mtick
from simplexTheory import *

filename = "Paper Figures/phasePlotIndividualContagionHealing"
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
ax.set_xlabel(r"$\beta_2$", fontsize=16)
ax.set_ylabel(r"$\beta_3$", fontsize=16)
xticks = np.linspace(xMin, xMax, 5)
yticks = np.linspace(yMin, yMax, 5)
plt.xticks(xticks, fontsize=12)
plt.yticks(yticks, fontsize=12)
ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))
ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))
plt.tight_layout()
plt.savefig('phasePlot.png',bbox_inches='tight',dpi = 300)
plt.show()
