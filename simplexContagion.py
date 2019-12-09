import random
import scipy
from scipy.sparse import csr_matrix
import numpy as np
import multiprocessing as mp


def microscopicSimplexSISDynamics(A, simplexList, simplexIndices, gamma, beta, alpha, x0, timesteps, dt):
    """Dynamical system"""
    n = np.shape(x0)[0]
    averageX = np.empty(timesteps)
    X = x0.copy()
    for i in range(timesteps-1):
        averageX[i] = np.mean(X)
        if averageX[i] == 0:
            X[random.randrange(n)] = 1

        # u = np.random.uniform(size=n)
        # v = np.random.uniform(size=n)
        pInfectByNeighbors = 1-np.power(1-beta*dt, A*X)
        for j in range(n):
            if X[j]==1:
                # heal
                X[j] = random.random() > gamma*dt
            else:
                if random.random() <= pInfectByNeighbors[j]:
                    # infect by neighbor
                    X[j] = 1
                else:
                    # don't run if there is zero contribution from simplices
                    if alpha > 0:
                        # infect by simplex
                        numInfectedSimplices = 0
                        for index in simplexIndices[j]:
                            if sum([X[member] for member in simplexList[index]]) >= 2:
                                numInfectedSimplices = numInfectedSimplices + 1
                        pInfectBySimplex = 1-(1-alpha*dt)**numInfectedSimplices
                        if random.random() <= pInfectBySimplex:
                            X[j] = 1

    averageX[-1] = np.mean(X)
    return averageX, X


def generateSISEquilibria(A, simplexList, simplexIndices, gamma, beta, alpha, x0, timesteps, dt, avgLength, verbose=True):
    equilibria = np.empty([len(alpha), len(beta)])
    for i in range(len(alpha)):
        x = x0.copy()
        for j in range(len(beta)):
            [sol, x] = microscopicSimplexSISDynamics(A, simplexList, simplexIndices, gamma, beta[j], alpha[i], x, timesteps, dt)
            equilibria[i,j] = np.mean(sol[-avgLength:-1])
            if verbose:
                print('alpha='+str(alpha[i])+', beta='+str(beta[j]), flush=True)
    return equilibria


def generateSISEquilibriaParallelized(A, simplexList, simplexIndices, gamma, beta, alpha, x0, timesteps, dt, avgLength, numProcesses):
    argList = []
    for alphaVal in alpha:
        argList.append((A, simplexList, simplexIndices, gamma, beta, alphaVal, x0, timesteps, dt, avgLength))
    with mp.Pool(processes=numProcesses) as pool:
        equilibria = pool.starmap(runOneCurve, argList)

    return equilibria


def runOneCurve(A, simplexList, simplexIndices, gamma, beta, alpha, x, timesteps, dt, avgLength, verbose=True):
    equilibria = []
    for i in range(len(beta)):
        [sol, x] = microscopicSimplexSISDynamics(A, simplexList, simplexIndices, gamma, beta[i], alpha, x, timesteps, dt)
        equilibria.append(np.mean(sol[-avgLength:-1]))
        if verbose:
            print('alpha='+str(alpha)+', beta='+str(beta[i]), flush=True)
    return equilibria


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
