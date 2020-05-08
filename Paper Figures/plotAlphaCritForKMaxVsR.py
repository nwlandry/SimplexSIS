import sys
sys.path.insert(1, 'C:/Users/nicho/Documents/GitHub/SimplexSIS')

import pickle
import matplotlib.pyplot as plt
import numpy as np

filename = "Paper Figures/alphaCritGridUncorrelated"
filename = "alphaCrit05072020-230602"
with open(filename, 'rb') as file:
    data = pickle.load(file)

xMin = data[0]
xMax = data[1]
yMin = data[2]
yMax = data[3]
alphaCritGrid = data[4]
firstOrderAlphaCritGrid = data[5]

# plt.figure()
# c = plt.imshow(np.flipud(alphaCritGrid), interpolation="spline16", cmap="Reds", extent=[xMin, xMax, yMin, yMax], aspect="auto")
# cbar = plt.colorbar(c)
# cbar.set_label(r"$\alpha_{crit}$", rotation=90)
# plt.xlabel("Power-Law Exponent")
# plt.ylabel("Maximum degree")
# plt.plot()
# plt.show()

plt.figure()
x = np.linspace(xMin, xMax, np.size(alphaCritGrid, axis=1))
y = np.linspace(yMin, yMax, np.size(alphaCritGrid, axis=0))
c = plt.contour(x, y, alphaCritGrid, 6)
plt.clabel(c, inline=1, fontsize=10)
plt.xlabel("Power-Law Exponent")
plt.ylabel("Maximum degree")
plt.plot()
plt.show()

plt.figure()
x = np.linspace(xMin, xMax, np.size(alphaCritGrid, axis=1))
y = np.linspace(yMin, yMax, np.size(alphaCritGrid, axis=0))
c = plt.contour(x, y, firstOrderAlphaCritGrid-alphaCritGrid, 5)
plt.clabel(c, inline=1, fontsize=10)
plt.xlabel("Power-Law Exponent")
plt.ylabel("Maximum degree")
plt.plot()
plt.show()
