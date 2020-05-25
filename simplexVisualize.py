import matplotlib.pyplot as plt
import numpy as np
import simplexTheory

def plotTheoreticalInfectionCurves(gamma, beta, alphaCritFractions, alphaCrit, degreeHist, meanSimplexDegree=None, isDegreeCorrelated=True, digits=4):
    if meanSimplexDegree == None:
        meanSimplexDegree = sum([k*prob for k, prob in degreeHist])
    plt.figure()
    color = ["green","blue","red","black","magenta"]
    i = 0
    for frac in alphaCritFractions:
        pointNum = 0
        for betaVal in beta:
            roots = simplexTheory.solveEquilibrium(gamma, betaVal, frac*alphaCrit, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits)
            for root in roots:
                if pointNum == 0:
                    plt.scatter(betaVal, root, color=color[i], label=str(alphaCritFractions[i])+r"$\beta_3^{(crit)}$")
                else:
                    plt.scatter(betaVal, root, color=color[i])
                pointNum = pointNum + 1
        i = i + 1
    plt.legend()
    plt.xlabel(r'$\beta_2$')
    plt.ylabel('infected average')
    plt.show()

def plotTheoreticalAndSimInfectionCurves(equilibrium, gamma, beta, alpha, degreeHist, meanSimplexDegree=None, isDegreeCorrelated=True, numTheoryPoints=50, digits=4):
    plt.figure()
    betaTheory = np.linspace(min(beta), max(beta), numTheoryPoints)
    for betaVal in betaTheory:
        roots = simplexTheory.solveEquilibrium(gamma, betaVal, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits)
        print(roots)
        for root in roots:
            plt.scatter(betaVal, root, color='black')

    plt.plot(beta, equilibrium, 'o-', color='blue')
    plt.xlabel(r'$\beta_2$')
    plt.ylabel('infected average')
    plt.show()

def calculateBistability(equilibrium, beta, option='area'):
    # These depend on there being an odd number of points
    if option == 'area':
        sumLower = 0
        sumUpper = 0
        for i in range(int((len(equilibrium)+1)/2)):
            sumLower = sumLower + (equilibrium[i+1]+equilibrium[i])/2.0*(beta[i+1] - beta[i])
            sumUpper = sumUpper + (equilibrium[-i-2]+equilibrium[-i-1])/2.0*(beta[-i-2] - beta[-i-1])
        return sumUpper - sumLower

    # Calculates the infinity norm between up and down curves
    elif option == 'infinity':
        return max([abs(equilibrium[i]-equilibrium[-i-1]) for i in range(int((len(equilibrium)+1)/2))])
    else:
        print("Please select a valid option")
