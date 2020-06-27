import sys
sys.path.insert(1, 'C:/Users/nicho/Documents/GitHub/SimplexSIS')

import pickle
import matplotlib.pyplot as plt
import numpy as np

filename = "Paper Figures/criticalRatioDegreeCorrelated"
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
plt.xlabel(r"$r$", fontsize=16)
plt.ylabel(r"$k_{max}$", fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig('criticalRatio.png',bbox_inches='tight',dpi = 300)
plt.show()

plt.figure()
x = np.linspace(xMin, xMax, np.size(alphaCritGrid, axis=1))
y = np.linspace(yMin, yMax, np.size(alphaCritGrid, axis=0))
c = plt.contour(x, y, np.divide(firstOrderAlphaCritGrid-alphaCritGrid,alphaCritGrid), 10)
plt.clabel(c, inline=1, fontsize=10)
plt.xlabel(r"$r$", fontsize=16)
plt.ylabel(r"$k_{max}$", fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig('criticalRatioError.png',bbox_inches='tight',dpi = 300)
plt.show()
