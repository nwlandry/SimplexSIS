from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.polynomial import polyroots

gamma = 2

#betaList = [0.0249]
alpha = 0.075
avgK = 75
avgSquaredK = 5833.3
avgCubedK = 468750
avgFourthK = 38750000
avgFifthK = 3281250000
avgSixthK = 1984375000000/7.0
numSteps = 0
step = 0.01
betaCrit = avgK/avgSquaredK*gamma
betaList = np.linspace(0, 0.1, 1000)
print(betaCrit)

def solveEquilbrium(gamma, betaList, alpha, avgK, avgSquaredK, avgCubedK, avgFourthK, avgFifthK, avgSixthK):
    equilibrium = list()
    for beta in betaList:
        #equilibrium.append(fsolve(equilibriumFunction, np.linspace(0,5,50),  args=(gamma, beta, alpha, avgK, avgSquaredK, avgCubedK, avgFourthK, avgFifthK)))
        roots = polyroots(equilibriumMatrix(gamma, beta, alpha, avgK, avgSquaredK, avgCubedK, avgFourthK, avgFifthK, avgSixthK))
        realRoots = list()
        for root in roots:
            if not np.iscomplex(root):
                realRoots.append(root)
        if len(realRoots) != 0:
            equilibrium.append(realRoots)
    return equilibrium

def equilibriumFunction(V, gamma, beta, alpha, avgK, avgSquaredK, avgCubedK, avgFourthK, avgFifthK):
    constant = beta/gamma*avgSquaredK/avgK - 1
    linearTerm = np.multiply(alpha/gamma*avgSquaredK/avgK - beta**2/gamma**2*avgCubedK/avgK,V)
    quadraticTerm = np.multiply(-2*alpha*beta/gamma**2*avgCubedK/avgK + beta**3/gamma**3*avgFourthK/avgK,np.power(V,2))
    cubicTerm = np.multiply(alpha**2/gamma**2*avgCubedK/avgK - 3*alpha*beta**2/gamma**3*avgFourthK/avgK + beta**4/gamma**4*avgFifthK/avgK,np.power(V,3))
    # print(str(constant) + ' + ' + str(linearTerm) + 'x + ' + str(quadraticTerm) + 'x^2 + ' + str(cubicTerm) + 'x^3')
    return constant + linearTerm + quadraticTerm + cubicTerm

def equilibriumMatrix(gamma, beta, alpha, avgK, avgSquaredK, avgCubedK, avgFourthK, avgFifthK, avgSixthK):
    constantTerm = beta/gamma*avgSquaredK/avgK - 1
    linearTerm = alpha/gamma*avgSquaredK/avgK - beta**2/gamma**2*avgCubedK/avgK
    quadraticTerm = -2*alpha*beta/gamma**2*avgCubedK/avgK + beta**3/gamma**3*avgFourthK/avgK
    cubicTerm = alpha**2/gamma**2*avgCubedK/avgK - 3*alpha*beta**2/gamma**3*avgFourthK/avgK + beta**4/gamma**4*avgFifthK/avgK
    quarticTerm = 3*alpha**2*beta/gamma**3*avgFourthK/avgK - 4*alpha*beta**3/gamma**4*avgFifthK/avgK + beta**5/gamma**5*avgSixthK/avgK
    #quinticTerm = alpha**2/gamma**2*avgCubedK/avgK - 3*alpha*beta**2/gamma**3*avgFourthK/avgK + beta**4/gamma**4*avgFifthK/avgK

    # print(str(constantTerm) + ' + ' + str(linearTerm) + 'x + ' + str(quadraticTerm) + 'x^2 + ' + str(cubicTerm) + 'x^3')
    coefficientMatrix = np.zeros(5)
    coefficientMatrix[0] = constantTerm
    coefficientMatrix[1] = linearTerm
    coefficientMatrix[2] = quadraticTerm
    coefficientMatrix[3] = cubicTerm
    coefficientMatrix[4] = quarticTerm
    #coefficientMatrix[5]
    return coefficientMatrix

equilibrium = solveEquilbrium(gamma, betaList, alpha, avgK, avgSquaredK, avgCubedK, avgFourthK, avgFifthK, avgSixthK)
print(equilibrium)

plt.figure()
for i in range(len(equilibrium)):
    for root in equilibrium[i]:
        plt.scatter(betaList[i], root, s=2)
#plt.plot(betaList, equilibrium)
plt.ylim([-2, 2])
plt.xlabel(r'$\beta$')
plt.ylabel('V')
plt.show()
