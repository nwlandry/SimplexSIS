import sys
sys.path.insert(1, 'C:/Users/nicho/Documents/GitHub/SimplexSIS')

import pickle
import matplotlib.pyplot as plt
import numpy as np

filename = "Paper Figures/ratioUncorrelated"
with open(filename, 'rb') as file:
    data = pickle.load(file)

xMin = data[0]
xMax = data[1]
yMin = data[2]
yMax = data[3]
alphaCritGrid = data[4]
firstOrderAlphaCritGrid = data[5]

plt.figure()
x = np.linspace(xMin, xMax, np.size(alphaCritGrid, axis=1))
y = np.linspace(yMin, yMax, np.size(alphaCritGrid, axis=0))
c = plt.contour(x, y, alphaCritGrid, 10)
plt.clabel(c, inline=1, fontsize=10)
plt.xlabel(r"$r$", fontsize=18)
plt.ylabel(r"$k_{max}$", fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.show()

plt.figure()
x = np.linspace(xMin, xMax, np.size(alphaCritGrid, axis=1))
y = np.linspace(yMin, yMax, np.size(alphaCritGrid, axis=0))
c = plt.contour(x, y, np.divide(firstOrderAlphaCritGrid-alphaCritGrid,alphaCritGrid), 10)
plt.clabel(c, inline=1, fontsize=10)
plt.xlabel(r"$r$", fontsize=18)
plt.ylabel(r"$k_{max}$", fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.show()

#
# plt.figure()
# c = plt.imshow(np.flipud(alphaCritGrid), interpolation="None", cmap="Reds", extent=[xMin, xMax, yMin, yMax], aspect="auto")
# cbar = plt.colorbar(c)
# cbar.set_label(r"$\alpha_{crit}$", rotation=90)
# plt.xlabel("Power-Law Exponent")
# plt.ylabel("Maximum degree")
# plt.show()
#
#
# plt.figure()
# c = plt.imshow(np.flipud(np.divide(firstOrderAlphaCritGrid-alphaCritGrid,alphaCritGrid)), interpolation="None", cmap="Reds", extent=[xMin, xMax, yMin, yMax], aspect="auto")
# cbar = plt.colorbar(c)
# cbar.set_label(r"$\alpha_{crit}$", rotation=90)
# plt.xlabel("Power-Law Exponent")
# plt.ylabel("Maximum degree")
# plt.show()
