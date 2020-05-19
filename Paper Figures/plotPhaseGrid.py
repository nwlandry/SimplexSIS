import sys
sys.path.insert(1, 'C:/Users/nicho/Documents/GitHub/SimplexSIS')

import pickle
import matplotlib.pyplot as plt
import numpy as np

filename = "Paper Figures/phasePlotIndividualDependent"
with open(filename, 'rb') as file:
    data = pickle.load(file)

xMin = data[0]
xMax = data[1]
yMin = data[2]
yMax = data[3]
phaseGrid = data[4]

plt.figure()
c = plt.imshow(np.flipud(phaseGrid), interpolation="none", cmap="Reds", extent=[xMin, xMax, yMin, yMax], aspect="auto", vmin=1, vmax=3)
plt.xlabel(r"$\beta_2$", fontsize=18)
plt.ylabel(r"$\beta_3$", fontsize=18)
plt.plot()
plt.show()
