
This file outlines how to reproduce the results from the real data analysis (section 3.3 and appendix A7). 

# Step 1: Construct knockoffs 

Knockoffs can be constructed following [Sesia et al. (2021)](https://www.pnas.org/doi/10.1073/pnas.2105841118). Details on the knockoff construction can be found [here](https://msesia.github.io/knockoffgwas/). The software for knockoff construction is available [here](https://github.com/msesia/knockoffgwas). 

# Step 2: Construct cloaked data 

To construct the cloaked data, two steps are required. 

a) Run functions in **batch_scripts/julia_functions_swapping_gof**.

This automatically submits separate .sh scripts to cloak the data (randomly swap original and knockoff genotypes) separately by chromosome. Internally, *batch_scripts/julia_functions_swapping_gof* uses the following file (you do not need to run it separately): 
   - *utils/swap_chr.R*.

b) Run the script **batch_scripts/swap_UKB_chrm.sh**. 

This merges the chromosome specific files back together, orders them as in the original file and to produce goodness-of-fit (GOF) plots. Internally, *utils/swap_UKB_chrm.sh* uses the following files (you do not need to run them separately):
  - *utils/merge_chromosomes.sh* 
  - *utils/make_FBM.R* 
  - *utils/order_chr_merged_file.R*
  - *utils/knockoffs_gof_swapping.R*

# Step 3: Local knockoff filter 

Run the script **batch_scripts/swap_UKB_whr_interactions.sh**. This single script will automatically conduct the steps a) to e) below (these steps do not need to be run individually): 

a) Pre-screening:
  - *utils/lasso_swapped.R*

b) Create interactions based on pre-screened variables:
  - *utils/create_interactions.R*: Creates interactions separately for each covariate. 
  - *utils/merge_interactions.sh*, *make_FBM.R*: Merge interactions to combine interactions from multiple covariates.

c) Run lasso with interactions, pick top interactions per group:
  - *utils/lasso_interactions.R*

d) Run lasso within each partition:
  - *utils/lasso_sep_each_environment.R*: Submit different batch scripts for each environment.
  - *utils/run_lasso_function_environment.R*: The batch scripts submitted above use this script to run the lasso within each partition.

e) Run the global-KF method, and also run the lasso directly after the pre-screening: 
  - *utils/lasso_real_data_after_prescreen.R*

# Step 4: Collect results

- After running *swap_UKB_wnb_interactions.sh*, run **utils/combine_lasso_environments.R** to collect importance statistics from each partition and pass through the knockoff to obtain aLKF rejections. 

- The script **utils/compare_rejections_aLKF_globalKF_fixedKF.R** compares the rejections between the aLKF, global KF and the Fixed-LKF for the Venn diagramm and table in appendix A7 and table in section 3.3.
- The script **utils/manhattan_plot.R** creates the Manhattan plot in section 3.3. 

- Note: To run the separate analysis for males and females, submit the script **utils/kf_male_female.sh**, which uses (*utils/male_female_whr.R*).

