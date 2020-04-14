import matplotlib.pyplot as plt
import numpy as np
import simplexTheory

def plotTheoreticalInfectionCurves(gamma, beta, alphaCritFractions, alphaCrit, degreeHist, meanSimplexDegree=None, isIndependent=False, digits=4):
    if meanSimplexDegree == None:
        meanSimplexDegree = sum([k*prob for k, prob in degreeHist])
    plt.figure()
    color = ["green","blue","red","black","magenta"]
    i = 0
    for frac in alphaCritFractions:
        pointNum = 0
        for betaVal in beta:
            roots = simplexTheory.solveEquilibrium(gamma, betaVal, frac*alphaCrit, degreeHist, meanSimplexDegree=meanSimplexDegree, isIndependent=isIndependent, digits=digits)
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

def plotTheoreticalAndSimInfectionCurves(equilibrium, gamma, beta, alpha, degreeHist, meanSimplexDegree=None, isIndependent=False, numTheoryPoints=50, digits=4):
    plt.figure()
    betaTheory = np.linspace(min(beta), max(beta), numTheoryPoints)
    for betaVal in betaTheory:
        roots = simplexTheory.solveEquilibrium(gamma, betaVal, alpha, minDegree, maxDegree, meanSimplexDegree, degreeHist=degreeHist, isIndependent=isIndependent, type=type, r=r, digits=digits)
        print(roots)
        for root in roots:
            plt.scatter(betaVal, root, color='black')

    plt.plot(beta, equilibrium, 'o-', color='blue')
    plt.xlabel(r'$\beta_2$')
    plt.ylabel('infected average')
    plt.show()

def calculateHysteresis(equilibrium, beta, option='area'):
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

# def findDiffMatrix(equilibria):
#     # These depend on there being an odd number of points
#     n = np.size(equilibria, axis=0)
#     m = int(0.5*(np.size(equilibria, axis=1)+1))
#     A = np.zeros([n, m])
#     for i in range(n):
#         for j in range(m):
#             A[i, j] = abs(equilibria[i][j]-equilibria[i][-j-1])
#     return A
#
# def findDiffMatrixAndGrid(equilibria, lambdaNetwork, lambdaSimplex):
#     # These depend on there being an odd number of points
#     n = np.size(equilibria, axis=0)
#     m = int(0.5*(np.size(equilibria, axis=1)+1))
#     A = np.zeros([n, m])
#     X = np.zeros([n, m])
#     Y = np.zeros([n, m])
#     for i in range(n):
#         for j in range(m):
#             A[i, j] = abs(equilibria[i][j]-equilibria[i][-j-1])
#             X[i, j] = lambdaNetwork[j]
#             Y[i, j] = lambdaSimplex[i]
#     return X, Y, A
#
#
# def findHysteresisBoundaries(equilibria, lambdaNetwork, lambdaSimplex, tolerance):
#     # These depend on there being an odd number of points
#     n = np.size(equilibria, axis=0)
#     m = int(0.5*(np.size(equilibria, axis=1)+1))
#     leftBoundary = []
#     rightBoundary = []
#     lambdaSimplexHys = []
#     for i in range(n):
#         left = max(lambdaNetwork)
#         right = 0
#         isHysteresis = False
#         for j in range(m):
#             if abs(equilibria[i][j]-equilibria[i][-j-1]) > tolerance:
#                 isHysteresis = True
#                 if lambdaNetwork[j] < left:
#                     left = lambdaNetwork[j]
#                 if lambdaNetwork[j] > right:
#                     right = lambdaNetwork[j]
#
#         if isHysteresis:
#             leftBoundary.append(left)
#             rightBoundary.append(right)
#             lambdaSimplexHys.append(lambdaSimplex[i])
#     return leftBoundary, rightBoundary, lambdaSimplexHys
#
# def plotTheoryAndSim(X, Y, A, leftBoundary, rightBoundary, lambdaSimplexHys):
#     plt.figure()
#     c = plt.pcolor(X, Y, A, cmap="Reds")
#     plt.plot(leftBoundary, lambdaSimplexHys, "k--")
#     plt.plot(rightBoundary, lambdaSimplexHys, "k--")
#     cbar = plt.colorbar(c)
#     cbar.set_label('Distance between fixed point solutions', rotation=90)
#     plt.xlabel(r"$\lambda$")
#     plt.ylabel(r"$\lambda_{\alpha}$")
#     plt.show()


# def plotSimHeatmapAndTheory(alpha, beta, A, leftBoundary, rightBoundary, alphaHys):
#     plt.figure()
#     xMin = min(beta)
#     xMax = max(beta)
#     yMin = min(alpha)
#     yMax = max(alpha)
#     c = plt.imshow(np.flipud(A), interpolation="spline16", cmap="Reds", extent=[xMin, xMax, yMin, yMax], aspect="auto")
#     plt.plot(leftBoundary, alphaHys, "k--")
#     plt.plot(rightBoundary, alphaHys, "k--")
#     cbar = plt.colorbar(c)
#     cbar.set_label('Distance between fixed point solutions', rotation=90)
#     plt.xlabel(r"$\beta$")
#     plt.ylabel(r"$\alpha$")
#     plt.show()
#
# def plotCompareSimHeatmapsAndTheory(alpha1, beta1, A1, leftBoundary1, rightBoundary1, alphaHys1, alpha2, beta2, A2, leftBoundary2, rightBoundary2, alphaHys2):
#     xMin1 = min(beta1)
#     xMax1 = max(beta1)
#     yMin1 = min(alpha1)
#     yMax1 = max(alpha1)
#     xMin2 = min(beta2)
#     xMax2 = max(beta2)
#     yMin2 = min(alpha2)
#     yMax2 = max(alpha2)
#
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#
#     c1 = ax.imshow(np.flipud(A1), interpolation="spline16", cmap="Reds", extent=[xMin1, xMax1, yMin1, yMax1], aspect="auto", alpha=0.4)
#     ax.plot(leftBoundary1, alphaHys1, "r--", label="Theory (Independent)")
#     ax.plot(rightBoundary1, alphaHys1, "r--")
#
#     c2 = ax.imshow(np.flipud(A2), interpolation="spline16", cmap="Blues", extent=[xMin2, xMax2, yMin2, yMax2], aspect="auto", alpha=0.4)
#     ax.plot(leftBoundary2, alphaHys2, "b--", label="Theory (Dependent)")
#     ax.plot(rightBoundary2, alphaHys2, "b--")
#
#     plt.xlim([max(xMin1,xMin2), min(xMax1, xMax2)])
#     plt.ylim([max(yMin1,yMin2), min(yMax1, yMax2)])
#     plt.xlabel(r"$\beta$")
#     plt.ylabel(r"$\alpha$")
#
#     ax.legend(loc="lower left")
#     # cbar1 = plt.colorbar(c1)
#     # cbar1.set_label('Simulation (Independent 1)', rotation=90)
#     # cbar1.solids.set_edgecolor("face")
#     # plt.draw()
#     # # Adding the colorbar
#     # cbaxes = fig.add_axes([0.01, 0.1, 0.03, 0.8])  # This is the position for the colorbar
#     # cbar2 = plt.colorbar(c2, cax = cbaxes, pad=0.1)
#     # cbar2.set_label('Simulation (Independent 2)')
#     # cbar2.solids.set_edgecolor("face")
#     # plt.draw()
#
#     plt.show()
#
# def plotCompareSimHeatmapsAndCritPoints(alpha1, beta1, A1, alphaCrit1, betaCrit1, alpha2, beta2, A2, alphaCrit2, betaCrit2):
#     xMin1 = min(beta1)
#     xMax1 = max(beta1)
#     yMin1 = min(alpha1)
#     yMax1 = max(alpha1)
#     xMin2 = min(beta2)
#     xMax2 = max(beta2)
#     yMin2 = min(alpha2)
#     yMax2 = max(alpha2)
#
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#
#     c1 = ax.imshow(np.flipud(A1), interpolation="spline16", cmap="Reds", extent=[xMin1, xMax1, yMin1, yMax1], aspect="auto", alpha=0.4)
#     ax.scatter(betaCrit1, alphaCrit1, label="Theory (Independent)", color="red")
#
#     c2 = ax.imshow(np.flipud(A2), interpolation="spline16", cmap="Blues", extent=[xMin2, xMax2, yMin2, yMax2], aspect="auto", alpha=0.4)
#     ax.scatter(betaCrit2, alphaCrit2, label="Theory (Dependent)", color="blue")
#
#     plt.xlim([max(xMin1,xMin2), min(xMax1, xMax2)])
#     plt.ylim([max(yMin1,yMin2), min(yMax1, yMax2)])
#     plt.xlabel(r"$\beta$")
#     plt.ylabel(r"$\alpha$")
#
#     ax.legend(loc="lower left")
#     # cbar1 = plt.colorbar(c1)
#     # cbar1.set_label('Simulation (Independent 1)', rotation=90)
#     # cbar1.solids.set_edgecolor("face")
#     # plt.draw()
#     # # Adding the colorbar
#     # cbaxes = fig.add_axes([0.01, 0.1, 0.03, 0.8])  # This is the position for the colorbar
#     # cbar2 = plt.colorbar(c2, cax = cbaxes, pad=0.1)
#     # cbar2.set_label('Simulation (Independent 2)')
#     # cbar2.solids.set_edgecolor("face")
#     # plt.draw()
#
#     plt.show()
