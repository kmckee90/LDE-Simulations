# LDE-Simulations
Simulations and tools for latent differential equations in R
Author: Kevin McKee
First commit: 9/23/2017

The main folder is the 'workspace' of the project and contains all the final simulation scripts and plot compilations.

Demo.R is set up to generate and analyze one simulated series at a time and is good for demonstration purposes.
compileVis.R is where I generate pdfs and other presentations.
Files starting with RunSim are the scripts to run complete simulations for each model and algorithm.

How to use:
Dependencies must be loaded in order, as appears in Demo.R or the RunSim files:

1. libraries, source everything in /lib/, as well as the resources you want to use in NL or MS.
  lib/NL/N_* (N-linked) scripts allow you to generate data for, and model, any sized network of stochastic differential equations.
  lib/MS/MS_* (Multi-scale) scripts generate data for, and model, univariate stochastic differential equations with dynamic behavior across two timescales.
  lib/MS/MSD_* (Multi-scale Discontinuous) is a variation of the multi-scale model where points are averaged with their neighbors and data is generated and organized to allow regular discontinuities. e.g. Discontinuities separating days, but continuity between day averages.
2. Set up a list() of default values with the correct names as shown in the Demo or RunSim files. Options and variables depend on the model being tested. (e.g. multiscale uses 2 embedding dimensions, N-linked only uses 1)
3. Load in order: matrices -> model -> constraints. 
3. For simulations, set up the matrices of parameters, means, and initial conditions, as found in the RunSim files.
4. Generate data using the corresponding genData file N_genData.R, MS_genData.R, MSD_genData.R
4b. simPar files simply take the setting of interest, the range of values to simulate for that setting, the matrices of true values and initial conditions, and your default values. They output a matrix of the final estimates, their standard errors, and the true values used to generate the data.
5. Time-embed data according to the specifications of the model with gllaembed().
6. Enter your embedded data into the model with model$data<-mxData(data, type="raw")
7. Run the model with with mxRun() or eventAnalysis(). 

eventAnalysis() outputs a list with:
$fit: the final fit of the OpenMx model
$exclusions: a vector or list of vectors of row indices from which to exclude data in the final iteration of model fitting.
$info: a matrix with columns representing the marginal -2log likelihoods of the data points of each series.

fitODE() takes an output from eventAnalysis, uses the model and exclusion indexes to calculate a purely deterministic ODE of best fit, which can then be plotted over the data. fitODE() outputs the mean squared error, mean absolute error, the predicted trajectory, and the estimated means, and is called by the option 'prediction=T' in eventAnalysis.

expectedValsML.R (currently unused) generates ML factor scores for the latent variables of the LDE. Good for examining how the LDE has interpreted process noise and measurement noise versus the deterministic components.

plotSim.R is for easily plotting simulation results. Will contain multiple plotting functions.
