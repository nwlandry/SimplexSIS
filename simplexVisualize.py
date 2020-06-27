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

def plotTheoreticalAndSimInfectionCurves(axisHandle, equilibrium, gamma, beta, alpha, degreeHist, meanSimplexDegree=None, isDegreeCorrelated=True, numTheoryPoints=50, digits=4):
    lineWidth = 2
    meanDegree = simplexTheory.computeMeanPowerOfDegreeFromHist(degreeHist, 1)
    meanSquaredDegree = simplexTheory.computeMeanPowerOfDegreeFromHist(degreeHist, 2)
    plt.figure()
    betaTheory = np.linspace(min(beta), max(beta), numTheoryPoints)
    theorySection1 = list() # zero stable
    theorySection2 = list() # zero unstable
    theorySection3 = list() # unstable branch
    theorySection4 = list() # stable epidemic
    ratioSection1 = list() # zero stable
    ratioSection2 = list() # zero unstable
    ratioSection3 = list() # unstable branch
    ratioSection4 = list() # stable epidemic

    for betaVal in betaTheory:
        roots = simplexTheory.solveEquilibrium(gamma, betaVal, alpha, degreeHist, meanSimplexDegree=meanSimplexDegree, isDegreeCorrelated=isDegreeCorrelated, digits=digits)
        if len(roots) == 1:
            theorySection1.append(roots[0])
            ratioSection1.append(betaVal/meanDegree*meanSquaredDegree/gamma)
        if len(roots) == 2:
            roots = sorted(roots)
            theorySection2.append(roots[0])
            ratioSection2.append(betaVal/meanDegree*meanSquaredDegree/gamma)
            theorySection4.append(roots[1])
            ratioSection4.append(betaVal/meanDegree*meanSquaredDegree/gamma)
        if len(roots) == 3:
            roots = sorted(roots)
            theorySection1.append(roots[0])
            ratioSection1.append(betaVal/meanDegree*meanSquaredDegree/gamma)
            theorySection3.append(roots[1])
            ratioSection3.append(betaVal/meanDegree*meanSquaredDegree/gamma)
            theorySection4.append(roots[2])
            ratioSection4.append(betaVal/meanDegree*meanSquaredDegree/gamma)
    axisHandle.plot(ratioSection1, theorySection1, 'k-',linewidth=lineWidth)
    axisHandle.plot(ratioSection2, theorySection2, 'k--',linewidth=lineWidth)
    axisHandle.plot(ratioSection3, theorySection3, 'k--',linewidth=lineWidth)
    axisHandle.plot(ratioSection4, theorySection4, 'k-',linewidth=lineWidth)
    if min(theorySection4) > 0.05 and len(theorySection3) <= 1:
        axisHandle.plot([ratioSection4[0], ratioSection4[0]], [0, theorySection4[0]], 'k--',linewidth=lineWidth)

    lineSegments = list()
    ratioSegments = list()
    currentLine = list()
    currentRatio = list()
    for i in range(len(equilibrium)-1):
        currentLine.append(equilibrium[i])
        currentRatio.append(beta[i]/meanDegree*meanSquaredDegree/gamma)
        if abs(equilibrium[i+1]-equilibrium[i])>0.1:
            if equilibrium[i] < equilibrium[i+1]:
                lineSegments.append(currentLine)
                ratioSegments.append(currentRatio)
                lineSegments.append([equilibrium[i], 2*equilibrium[i+1]-equilibrium[i+2]])
                ratioSegments.append([beta[i]/meanDegree*meanSquaredDegree/gamma,beta[i]/meanDegree*meanSquaredDegree/gamma])
                currentLine = [2*equilibrium[i+1]-equilibrium[i+2]]
                currentRatio = [beta[i]/meanDegree*meanSquaredDegree/gamma]
            else:
                lineSegments.append(currentLine)
                ratioSegments.append(currentRatio)
                lineSegments.append([equilibrium[i], equilibrium[i+1]])
                ratioSegments.append([beta[i]/meanDegree*meanSquaredDegree/gamma,beta[i]/meanDegree*meanSquaredDegree/gamma])
                currentLine = [equilibrium[i+1]]
                currentRatio = [beta[i]/meanDegree*meanSquaredDegree/gamma]
    lineSegments.append(currentLine)
    ratioSegments.append(currentRatio)

    for i in range(len(lineSegments)):
        if len(lineSegments[i]) == 2:
            axisHandle.plot(ratioSegments[i], lineSegments[i], '--', color='blue',linewidth=lineWidth)
        else:
            axisHandle.plot(ratioSegments[i], lineSegments[i], 'o-', color='blue',linewidth=lineWidth)
    axisHandle.set_xlim([max(0,min(beta)/meanDegree*meanSquaredDegree/gamma-0.01), max(beta)/meanDegree*meanSquaredDegree/gamma+0.01])
    axisHandle.set_ylim([-0.005, 0.7])# max(theorySection4)+0.01])
    plt.tight_layout()
    plt.close()
    # plt.xlabel(r'$\beta_2/\beta_2^c$', fontsize=16)
    # plt.ylabel(r'$U$', fontsize=16)
    # plt.tight_layout()
    # plt.show()

def calculateBistability(equilibrium, beta, option='infinity'):
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
