
# UKB GWAS
#
# Class: R script
#
# This script applies the knockoff filter on the knockoff statistics from each partition to get rejections 
# of the aLKF 

library(tidyverse)

# set parameters
my_population = "british"
resolution = 4
pheno.name = "whr"
fdr = 0.1
covar_identifiyer = "sex" 
suffix = ""
num.int = 1

# path where knockoff statistics are stored
stats.path = "path where to store results and knockoff statistics"

# path where information on relevant loci within each partition are stored
list_for_combine_lasso = readRDS(paste0("lasso_each_env/combination_list_",
                                                my_population, "_res", resolution, "_", 
                                        pheno.name, "_cov_", covar_identifiyer, "_", suffix, "_", num.int))


unique_interaction_vars = list_for_combine_lasso$unique_interaction_vars
all_environments = list_for_combine_lasso$all_environments

# combine knockoff statistics
all_stats = c()
file_does_not_exist = c()
for(i in 1:length(unique_interaction_vars)) {
  
  unique_environments = all_environments[[i]]
  for(e in unique_environments) {
    filename = paste0(stats.path, "lasso_env_stats_",my_population, "_res", resolution, "_", pheno.name, 
                      "_cov_", covar_identifiyer, "_",
                      paste(unique_interaction_vars[[i]], collapse = "."), "_env", e, "_", suffix, "_", num.int,".txt")
    
    
    print(paste0("Does this file exist? ", file.exists(filename)))
    
    if(file.exists(filename)) {
      Stats_filtered_file <-  read_delim(filename)
      all_stats = rbind(all_stats, Stats_filtered_file)
    } else {
      tmp_file = data.frame(covar = NA, env = NA)
      tmp_file[1, 1] <- paste(unique_interaction_vars[[i]], collapse = ".")
      tmp_file[1, 2] <- e
      
      file_does_not_exist = rbind(file_does_not_exist, tmp_file)
    }
    
  }

}

# check that all files exist
print(file_does_not_exist)

knockoff.threshold <- function(W, fdr=0.10, offset=1) {
  if(offset>1 | offset<0) {
    stop('Input offset must be between 0 or 1')
  }
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}
knockoff.filter <- function(Stats, fdr=0.1, offset=1) {
  W.thres <- knockoff.threshold(Stats$W, fdr=fdr, offset=offset)
  Selected <- Stats %>% filter(W >= W.thres)
  return(Selected)
}

# knockoff filter 
selections <- all_stats %>% knockoff.filter(fdr=fdr, offset=1) %>% 
  arrange(desc(W)) 

print(summary(all_stats$W))
print(summary(selections$W))
print(sum(all_stats$W > 0))
dim(selections)

# merge in BP information by group and save selections and statistics
min_max_BP_by_group.file <- sprintf("%s_group_bp_min_max.txt", paste0(stats.path, my_population, 
                                                                      "_res", resolution, "_", pheno.name, "_cov_", covar_identifiyer))
min_max_BP_by_group  = read_delim(min_max_BP_by_group.file, delim=" ")

selections = selections %>% 
  left_join(min_max_BP_by_group, by = c("CHR", "Group"))

all_stats = all_stats %>% 
  left_join(min_max_BP_by_group, by = c("CHR", "Group"))

write.csv(selections, 
          paste0(stats.path, "lasso_selections_interactions_weight0.25_" ,my_population, "_res",
                 resolution, "_", pheno.name, "_cov_", covar_identifiyer, "_", suffix, "_", num.int ,".csv"), row.names = FALSE)

write.csv(all_stats, 
          paste0(stats.path, "lasso_statistics_interactions_weight0.25_",
                 my_population, "_res", resolution, "_", pheno.name, "_cov_", covar_identifiyer, "_", suffix, "_", num.int , ".csv"), row.names = FALSE)



# combine coefficients 
all_betas = c()
for(i in 1:length(unique_interaction_vars)) {

  print(i)
  
  unique_environments = all_environments[[i]]
  for(e in unique_environments) {

      print(e)
    
    
    filename <- paste0(stats.path, "lasso_env_stats_",my_population, "_res", resolution, "_", pheno.name, "_cov_", covar_identifiyer, "_",
                       paste(unique_interaction_vars[[i]], collapse = "."), "_env", e, "_", suffix, "_", num.int , "_beta_filtered" ,".txt")

    print(filename)
    print(paste0("Does this file exist? ", file.exists(filename)))

    if(file.exists(filename)) {
      Beta_filtered_file <-  read_delim(filename) %>% mutate(Interaction = paste(unique_interaction_vars[[i]], collapse = "."), 
                                                             Environment = e)
      all_betas = rbind(all_betas, Beta_filtered_file)
    } else {
      tmp_file = data.frame(covar = NA, env = NA)
      tmp_file[1, 1] <- paste(unique_interaction_vars[[i]], collapse = ".")
      tmp_file[1, 2] <- e

      file_does_not_exist = rbind(file_does_not_exist, tmp_file)
    }

  }

}

# save coefficients
all_betas_filtered_selections = all_betas %>% filter(CHR_Group %in% selections$CHR_Group)

write.csv(all_betas_filtered_selections, 
          paste0(stats.path, "betas_coef_selections_interactions_weight0.25_" ,my_population, "_res",
                 resolution, "_", pheno.name, "_cov_", covar_identifiyer, "_", suffix, "_", num.int ,".csv"), row.names = FALSE)


