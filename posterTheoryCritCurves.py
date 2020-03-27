import simplexTheory
import visualizeData
import simplexContagion
import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

powers = np.linspace(3.001,6,100)
depAlphaCrit = list()
indepAlphaCrit = list()
n = 1000
meanDegreeTarget = 100
initialGuess = 0.1
gamma = 2


for power in powers:
    minDeg = fsolve(simplexTheory.avgOfPowerLawEqn, initialGuess,  args=(n, power, meanDegreeTarget))
    minDeg = max(minDeg, 1)
    print(minDeg)
    meanDegree = simplexTheory.meanPowerOfPowerLaw(minDeg, n, power, 1)
    meanSquaredDegree = simplexTheory.meanPowerOfPowerLaw(minDeg, n, power, 2)
    meanCubedDegree = simplexTheory.meanPowerOfPowerLaw(minDeg, n, power, 3)
    indepAlphaCrit.append(meanSquaredDegree/(meanDegree**3)*gamma)
    depAlphaCrit.append((meanDegree**2)*meanCubedDegree/(meanSquaredDegree**3)*gamma)


minUniform = 10
maxUniform = 30

meanDegreeUniform = 0.5*(minUniform + maxUniform)
meanSquaredDegreeUniform = 1.0/(3*maxUniform-3*minUniform)*(maxUniform**3 - minUniform**3)
meanCubedDegreeUniform = 1.0/(4*maxUniform-4*minUniform)*(maxUniform**4 - minUniform**4)

critUniformIndep = meanSquaredDegreeUniform/(meanDegreeUniform**3)*gamma
critUniform = (meanDegreeUniform**2)*meanCubedDegreeUniform/(meanSquaredDegreeUniform**3)*gamma

plt.figure()
plt.plot(powers, depAlphaCrit, linewidth=2, label=r"Critical $\alpha$ for dependent simplices (power law)")
plt.plot(powers, indepAlphaCrit, linewidth=2, label=r"Critical $\alpha$ for independent simplices (power law)")
plt.plot([min(powers), max(powers)], [critUniform, critUniform], linewidth=2, label=r"Critical $\alpha$ for dependent simplices (uniform)")
plt.plot([min(powers), max(powers)], [critUniformIndep, critUniformIndep], linewidth=2, label=r"Critical $\alpha$ for independent simplices (uniform)")
plt.legend()
plt.xlabel("Exponent of power law")
plt.ylabel(r"$\alpha_{crit}$")
plt.show()
