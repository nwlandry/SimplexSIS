# Overview of the SimplexSIS repository
Data and code used to generate results for the paper "The effect of heterogeneity on hypergraph contagion models"

## Library files
* simplexContagion.py: Implements the functions required to run the microscopic Markov model on a hypergraph
  - microscopicSimplexSISDynamics(): Implementation of a microscopic Markov model including 3-hyperedges. Outputs a population average time-series.
  - generateSISEquilibriaParallelized(): Runs several different equilibrium curves for each specified value of $\beta_3$ in parallel. I.e. it runs runOneCurve() for each $\beta_3$ value.
  - runOneCurve(): Runs a single equilibrium curve.
  - The rest of this file was not used in the paper
* simplexUtilities.py: Implements functions used to generate the adjacency matrix and the lists of hyperedges.
  - criticalBeta3(): Computes $\beta_3^c$ from the network properties (Methods in Section IV)
  - invCDFPowerLaw() and generatePowerLawDegreeSequence(): generate a power-law degree sequence.
  - generateUniformDegreeSequence(): Generates a uniformly distributed degree sequence.
  - generateConfigModelAdjacency(): Generates an adjacency matrix from a degree sequence using the configuration model
  - generateConfigModelSimplexList(): Generates a list of hyperedges using the configuration model
  - generateUniformSimplexList(): Generates a list of hyperedges that a uniformly distributed on the nodes
  - The rest of this file is not used in the paper
* simplexVisualize.py: Implements different functions for visualizing the mean field curves and bistability index
* simplexTheory.py: Implements functions dealing with the mean-field model
  - getPhase(): Gets the number of roots solving the mean-field equations
  - solveEquilibrium(): Solves the mean-field equation for given degree distribution, $\gamma$, $\beta_2$, and $\beta_3$
  - generateTheoreticalDegreeHist(): generates a degree list with probabilities from a given distribution
  - calculateTheoreticalCriticalAlpha(): Finds $\beta_3^c$ using methods described in Appendix A.3
  - The rest of this file were auxiliary functions or not used.
## General functions
* runModelInParallel.py: Runs the microscopic model in parallel
* outputPhasePlotBetaAlpha.py: Generates the phase plots for the mean-field models.
* outputRatioForKMaxVsRInParallel: Generates the contour plots of $\beta_3^c/\beta_2^c$ for $k_{max}$ vs. $r$.
## Paper Figures Folder: Functions and data used to plot the figures in the paper
