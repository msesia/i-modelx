
# Reproducing the genetic simulation results 


## Step 1: Run simulations

All simulations can be run using the functions provided in **batch_scripts/julia_functions_sim**. This allows to reproduce all simulation results from section 3.2 and appendix A6. Each function in *batch_scripts/julia_functions_sim* runs one method (aLKF, global-KF, LKF-naive, LKF_split, Fixed-LKF, uw-aLKF - see section 3 for further details). The functions should be run in the [julia language](https://julialang.org/). However, these functions only write and submit .sh files. The actual analyses are conducted in R. 

The files used in these functions are: 

- *sim_globalKF.R*: Simulation for the *global-KF*. 
- *sim_aLKF.R*: Simulation for the *aLKF*. 
- *sim_uwaLKF.R*: Simulation for the *uw-aLKF*. 
- *sim_LKF_naive.R*: Simulation for the *LKF-naive*. 
- *sim_LKF_split.R*: Simulation for the *LKF-split*. 
- *sim_fixed_LKF_1cov.R*: Simulation for the *Fixed-LKF* for one covariate.
- *sim_fixed_LKF_2cov.R*: Simulation for the *Fixed-LKF* for two covariates.

## Step 2: Make plots

After running all simulations, the script **utils/simulation_plots.R** allows to create the plots in sections 3.2 and appendix A6.


