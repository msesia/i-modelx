
This file outlines the scripts used for the real data analysis and explains the process. 

# Construct cloaked data 

a) Run functions in *batch_scripts/julia_functions_swapping_gof* to submit separate .sh scripts to cloak the data (randomly swap original and knockoff genotypes). This is conducted separately by chromosome. The file used is:
  - *utils/swap_chr.R*

b) Run *utils/swap_UKB_chrm.sh* to merge chromosome specific files back together, to order them as in the original file and to produce goodness-of-fit (GOF) plots. The files used in this script are:
  - *utils/merge_chromosomes.sh* 
  - *utils/make_FBM.R* 
  - *utils/order_chr_merged_file.R*
  - *utils/knockoffs_gof_swapping.R*

# Local knockoff filter #

Running the script *swap_UKB_wnb_interactions.sh* will run a) to e) below: 

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

After running *swap_UKB_wnb_interactions.sh*, run *utils/combine_lasso_environments.R* to collect importance statistics from each partition and pass through the knockoff to obtain aLKF rejections.

To run the separate analysis for males and females, submit the script *utils/kf_male_female.sh*, which uses (*utils/male_female_whr.R*).

