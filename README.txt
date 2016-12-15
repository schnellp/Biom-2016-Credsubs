FILES:
  - alzheimers.csv:   Contains cleaned data for the Alzheimer's
                      disease example.
  - bayes-linear.R:   Defines functions related to the methods
                      described in the paper.
  - decomp.R	      Produces estimate, standard error, and observation
    		      plots used in example analysis.
  - example.R:        Reproduces the Alzheimer's disease example analysis.
  - simulate.R:       Reproduces the simulation study.

TO RUN ALZHEIMER'S DISEASE ANALYSIS:
  - Packages required: MASS, Matrix, matrixcalc, matrixStats, mvtnorm,
    plotrix, vcd.
  - Make sure working directory is set to the directory
    containing the R scripts.
  - Run example.R.
  - To reproduce the decomposition plots, run decomp.R.

TO RUN SIMULATION STUDY:
  - Packages required: MASS, Matrix, matrixcalc, matrixStats, mvtnorm,
    BayesTree.
  - Make sure working directory is set to the directory
    containing the R scripts.
  - Run simulate.R. This may take several hours.
  - Simulation results are saved in tables-n.RData, where
    n is the sample size (default 40).

TO ANALYZE A NEW DATASET:
  - Load bayes-linear.R
  - Defaults: call blm() with syntax used by lm().
  - Custom prior: call blm(..., design=TRUE) to get the design
    matrix, and specify hyperparameters accordingly.
  - Create grid approximating the covariate space of interest.
  - Call the find.credible.subgroups function.
  - See example.R for an example use.
