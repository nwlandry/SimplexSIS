import matplotlib.pyplot as plt
import numpy as np
import simplexTheory

def plotTheoreticalInfectionCurves(gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, degreeSequence=None, isIndependent=False, type="power-law", r=4, digits=4):
    plt.figure()
    for betaVal in beta[0:int((len(beta)+1)/2)]:
        roots = simplexTheory.solveEquilibrium(gamma, betaVal, alpha, minDegree, maxDegree, meanSimplexDegree, degreeSequence=degreeSequence, isIndependent=isIndependent, type=type, r=r, digits=digits)
        for root in roots:
            plt.scatter(betaVal, root, color='black')

    # for betaVal in beta[0:int((len(beta)+1)/2)]:
    #     roots = simplexTheory.solveEquilibrium(gamma, betaVal, alpha, minDegree, maxDegree, meanSimplexDegree, degreeSequence=None, isIndependent=isIndependent, type=type, r=r, digits=digits)
    #     for root in roots:
    #         plt.scatter(betaVal, root, color='red')

    plt.xlabel(r'$\beta$')
    plt.ylabel('infected average')
    plt.show()

def plotTheoreticalAndSimInfectionCurves(equilibrium, gamma, beta, alpha, minDegree, maxDegree, meanSimplexDegree, degreeSequence=None, isIndependent=False, type="power-law", r=4, numTheoryPoints=50, digits=4):
    plt.figure()
    betaTheory = np.linspace(min(beta), max(beta), numTheoryPoints)
    for betaVal in betaTheory:
        roots = simplexTheory.solveEquilibrium(gamma, betaVal, alpha, minDegree, maxDegree, meanSimplexDegree, degreeSequence=degreeSequence, isIndependent=isIndependent, type=type, r=r, digits=digits)
        #roots = simplexTheory.solveEquilibriumOld(gamma, betaVal, alpha, minDegree, maxDegree, meanSimplexDegree, digits=4, isIndependent=isIndependent, type="power-law", r=4)
        for root in roots:
            plt.scatter(betaVal, root, color='black')

    for betaVal in beta[0:int((len(beta)+1)/2)]:
        roots = simplexTheory.solveEquilibrium(gamma, betaVal, alpha, minDegree, maxDegree, meanSimplexDegree, degreeSequence=None, isIndependent=isIndependent, type=type, r=r, digits=digits)
        for root in roots:
            plt.scatter(betaVal, root, color='red')

    plt.plot(beta, equilibrium, 'o-', color='blue')
    plt.xlabel(r'$\beta$')
    plt.ylabel('infected average')
    plt.show()

def plotSimHeatmapAndTheory(alpha, beta, A, leftBoundary, rightBoundary, alphaHys):
    plt.figure()
    xMin = min(beta)
    xMax = max(beta)
    yMin = min(alpha)
    yMax = max(alpha)
    c = plt.imshow(np.flipud(A), interpolation="spline16", cmap="Reds", extent=[xMin, xMax, yMin, yMax], aspect="auto")
    plt.plot(leftBoundary, alphaHys, "k--")
    plt.plot(rightBoundary, alphaHys, "k--")
    cbar = plt.colorbar(c)
    cbar.set_label('Distance between fixed point solutions', rotation=90)
    plt.xlabel(r"$\beta$")
    plt.ylabel(r"$\alpha$")
    plt.show()

def plotCompareSimHeatmapsAndTheory(alpha1, beta1, A1, leftBoundary1, rightBoundary1, alphaHys1, alpha2, beta2, A2, leftBoundary2, rightBoundary2, alphaHys2):
    xMin1 = min(beta1)
    xMax1 = max(beta1)
    yMin1 = min(alpha1)
    yMax1 = max(alpha1)
    xMin2 = min(beta2)
    xMax2 = max(beta2)
    yMin2 = min(alpha2)
    yMax2 = max(alpha2)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    c1 = ax.imshow(np.flipud(A1), interpolation="spline16", cmap="Reds", extent=[xMin1, xMax1, yMin1, yMax1], aspect="auto", alpha=0.4)
    ax.plot(leftBoundary1, alphaHys1, "r--", label="Theory (Independent)")
    ax.plot(rightBoundary1, alphaHys1, "r--")

    c2 = ax.imshow(np.flipud(A2), interpolation="spline16", cmap="Blues", extent=[xMin2, xMax2, yMin2, yMax2], aspect="auto", alpha=0.4)
    ax.plot(leftBoundary2, alphaHys2, "b--", label="Theory (Dependent)")
    ax.plot(rightBoundary2, alphaHys2, "b--")

    plt.xlim([max(xMin1,xMin2), min(xMax1, xMax2)])
    plt.ylim([max(yMin1,yMin2), min(yMax1, yMax2)])
    plt.xlabel(r"$\beta$")
    plt.ylabel(r"$\alpha$")

    ax.legend(loc="lower left")
    # cbar1 = plt.colorbar(c1)
    # cbar1.set_label('Simulation (Independent 1)', rotation=90)
    # cbar1.solids.set_edgecolor("face")
    # plt.draw()
    # # Adding the colorbar
    # cbaxes = fig.add_axes([0.01, 0.1, 0.03, 0.8])  # This is the position for the colorbar
    # cbar2 = plt.colorbar(c2, cax = cbaxes, pad=0.1)
    # cbar2.set_label('Simulation (Independent 2)')
    # cbar2.solids.set_edgecolor("face")
    # plt.draw()

    plt.show()

def plotCompareSimHeatmapsAndCritPoints(alpha1, beta1, A1, alphaCrit1, betaCrit1, alpha2, beta2, A2, alphaCrit2, betaCrit2):
    xMin1 = min(beta1)
    xMax1 = max(beta1)
    yMin1 = min(alpha1)
    yMax1 = max(alpha1)
    xMin2 = min(beta2)
    xMax2 = max(beta2)
    yMin2 = min(alpha2)
    yMax2 = max(alpha2)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    c1 = ax.imshow(np.flipud(A1), interpolation="spline16", cmap="Reds", extent=[xMin1, xMax1, yMin1, yMax1], aspect="auto", alpha=0.4)
    ax.scatter(betaCrit1, alphaCrit1, label="Theory (Independent)", color="red")

    c2 = ax.imshow(np.flipud(A2), interpolation="spline16", cmap="Blues", extent=[xMin2, xMax2, yMin2, yMax2], aspect="auto", alpha=0.4)
    ax.scatter(betaCrit2, alphaCrit2, label="Theory (Dependent)", color="blue")

    plt.xlim([max(xMin1,xMin2), min(xMax1, xMax2)])
    plt.ylim([max(yMin1,yMin2), min(yMax1, yMax2)])
    plt.xlabel(r"$\beta$")
    plt.ylabel(r"$\alpha$")

    ax.legend(loc="lower left")
    # cbar1 = plt.colorbar(c1)
    # cbar1.set_label('Simulation (Independent 1)', rotation=90)
    # cbar1.solids.set_edgecolor("face")
    # plt.draw()
    # # Adding the colorbar
    # cbaxes = fig.add_axes([0.01, 0.1, 0.03, 0.8])  # This is the position for the colorbar
    # cbar2 = plt.colorbar(c2, cax = cbaxes, pad=0.1)
    # cbar2.set_label('Simulation (Independent 2)')
    # cbar2.solids.set_edgecolor("face")
    # plt.draw()

    plt.show()
