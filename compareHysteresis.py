import pickle
import matplotlib.pyplot as plt
import numpy as np
import simplexUtilities
import simplexContagion


filenameUniform = 'equilibriaData10252019-225011'
filenamePowerLaw = 'equilibriaData11072019-102634'
with open(filenameUniform, 'rb') as file:
    data = pickle.load(file)
alphaUnif = data[0]
betaUnif = data[1]
equilibriaUnif = data[2]

with open(filenamePowerLaw, 'rb') as file:
    data = pickle.load(file)
alphaPL = data[0]
betaPL = data[1]
equilibriaPL = data[2]

hysteresisUnif = np.zeros([len(alphaUnif)])

for i in range(len(alphaUnif)):
    hysteresisUnif[i] = simplexContagion.calculateHysteresis(equilibriaUnif[i], betaUnif, option='area')

hysteresisPL = np.zeros([len(alphaPL)])
for i in range(len(alphaPL)):
    hysteresisPL[i] = simplexContagion.calculateHysteresis(equilibriaPL[i], betaPL, option='area')

plt.figure()
plt.plot(alphaUnif, hysteresisUnif, 'o-', label='Uniformly Distributed Degree')
plt.plot(alphaPL, hysteresisPL, 'o-', label='Power Law Degree Distribution')
plt.xlabel(r"$\alpha$")
plt.ylabel("Hysteresis")
plt.legend()
plt.show()
