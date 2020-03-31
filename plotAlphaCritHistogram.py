import simplexTheory
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from simplexTheory import *


filename = "alphaCritList03292020-212018"
with open(filename, 'rb') as file:
    data = pickle.load(file)

print(len(data))
plt.figure()
plt.hist(data, bins=40)
plt.xlabel(r"$\alpha_{crit}$")
plt.ylabel("Instances")
plt.plot()
plt.show()
