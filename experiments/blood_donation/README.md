# Reproducing the semi-synthetic blood donation experiment results

This folder contains the code needed to reproduce the results from *Appendix A4.2* using data from the blood donation randomized experiment.  

This folder does **not** contain any data files. We obtained access to anonymized data from the randomized study of blood donation incentives from:  

> Sun, T., G. Gao, and G. Z. Jin (2019). *Mobile messaging for offline group formation in prosocial activities: A large field experiment.* Management Science 65(6), 2717–2736.  

Interested parties can request access to these data by contacting the authors of that paper.  
The analysis was conducted in **R/4.2.0** and **Python/3.7.6**.

---

## Step 1: Construct knockoffs for treatments

Knockoffs for the binary treatment assignments are constructed using the Metropolized knockoff sampler of [Bates et al. (2020)](https://arxiv.org/abs/2002.08943).  

- Run the script **setup/submit_knockoffs.sh**.  
  This submits jobs (or runs locally if uncommented) to generate knockoffs for the 5 treatment variables, repeating the process across 100 random seeds.  

Internally, the following files are used (do not run them separately):  
- *setup/generate_knockoffs.py* — Python code implementing the Metropolized sampler.  
- *setup/generate_knockoffs.sh* — wrapper to call the Python script with a given seed.  

Outputs are stored in:  
- ../data/knockoffs/dt_donor_mar042018_simple_N5_seed${SEED}.csv


---

## Step 2: Simulate semi-synthetic outcomes

Because the real donation outcome is extremely imbalanced, we simulate semi-synthetic outcomes according to the causal logistic model in Appendix A4.2 (Eq. (A.4.2.1)).  

- Run the script **experiments/submit_generate_donations.sh**.  
  This calls **experiments/generate_donations.sh**, which in turn runs **experiments/generate_donations.R**.  

This script:  
- Reads in the real covariates and treatment + knockoff data.  
- Simulates semi-synthetic outcomes using the specified logistic model with parameters \(a=0.4\), \(b=0.2\).  
- Saves the generated outcomes into the data directory for use in subsequent experiments.  

---

## Step 3: Run knockoff filter experiments

The main numerical experiments are carried out with the adaptive local knockoff filter, along with benchmark methods.  

- Submit the script **experiments/submit_experiment.sh**.  
  This runs **experiments/experiment.sh**, which calls **experiments/experiment.R**.  

This script:  
- Loads the knockoff treatments and simulated outcomes.  
- Applies the adaptive local knockoff filter across multiple sample sizes.  
- Runs benchmark methods (global knockoff filter, naive local filter, data splitting).  
- Repeats the analysis across multiple seeds for independent replicates.  

Optional: **experiments/submit_experiment_2.sh** provides a variant using **experiment_2.R** with alternative configurations.  

---

## Step 4: Collect results and generate plots

- Run **experiments/make_plots.R** to compile the results across experiments and generate figures.  
  This script produces Figure A.4.2.1 (performance comparison of methods) and diagnostic plots for the knockoffs.  

Internally, the scripts *utils_ci.R* and *utils_ci_sharp.R* provide helper functions for confidence intervals and performance summaries.  

Outputs include:  
- Performance summaries (FDR, power, homogeneity).  
- Plots saved into the `figures/` directory, including **experiment_blood.pdf** (as in Appendix A4.2).  

---
