#!/usr/bin/env Rscript

# UKB GWAS
#
# Class: R script
#
# This script submits separate lasso jobs for each partition


# Input arguments
args <- commandArgs(trailingOnly=TRUE)
my_population       <- as.character(args[1])
resolution       <- as.integer(args[2])
pheno.name       <- as.character(args[3])
time_lasso <- as.integer(args[4])
memory_lasso <- as.integer(args[5])
dfmax <- as.integer(args[6])
ncores <- as.integer(args[7])
partition_sherlock <- as.character(args[8])
suffix <- as.character(args[9])
num.int <- as.numeric(args[10])
pheno.file <- as.character(args[11])
real_data_file <- as.character(args[12])

library(bigsnpr)
library(tidyverse)


set.seed(2020)

################## Define paths ################


covariate_file = paste0("covariates/covariates_", pheno.name, "_",
                        my_population ,"_sample_all_groups")


covariates_interactions = read.csv(covariate_file)

# get last part after "_"
get_last_string <- function(string) {
  parts <- strsplit(string, "_")[[1]]
  return(tail(parts, 1))
}
filtered_col_names <- colnames(covariates_interactions)[!grepl("eid", colnames(covariates_interactions)) & !grepl("^rev", colnames(covariates_interactions))]
covar_names_pretty <- sapply(filtered_col_names, get_last_string)
covar_identifiyer = paste0(sort(covar_names_pretty), collapse = "_")


stats.path = paste0("UKB_stats/", my_population, "_res", resolution, "_", pheno.name, "_cov_", covar_identifiyer, "_numint", num.int)

# load partition from lasso with interactions
partition.file <- sprintf("%s_lasso_interaction_partition.txt", stats.path)
partition = read_delim(partition.file, delim=" ") %>% mutate(CHR_Group = paste0(CHR, "_", Group))
unique_chr_group = unique(partition$CHR_Group)

print(paste0("rows in partition: ", nrow(partition)))
print(partition %>% group_by(Interaction) %>% summarize(count = n()))

partition_covariates <- function(Z, vars.split) {
  bitsToInt <- function(x) {
    packBits(rev(c(rep(FALSE, 32-length(x)%%32), as.logical(x))), "integer")
  }
  n = nrow(Z)
  ## Number of partitions
  K = length(vars.split)
  out = rep(-1, n)
  for(i in 1:n) {
    out[i] = 1 + bitsToInt(Z[i,vars.split])
  }
  return(out)
}

# calculate unique environments
all_interaction_vars = list()
for(r in unique_chr_group) {
  # filter partition to the group; arrange ensures that the ordering for all is the same
  chr_group_interaction = partition %>% dplyr::filter(CHR_Group == r) %>% arrange(Interaction)
  all_interaction_vars = c(all_interaction_vars, list(chr_group_interaction$Interaction))
}

unique_interaction_vars = unique(all_interaction_vars)
all_environments = lapply(unique_interaction_vars, function(x) partition_covariates(covariates_interactions %>% dplyr::select(-eid), x))

# find the chromosome group associated with each interaction
find_chr_group_by_interaction <- function(df, int) {

  given_interaction <- paste(sort(int), collapse = ", ")

  # Filter partition for matching interactions
  df %>%
    group_by(CHR_Group) %>%
    summarize(interactions = paste(sort(Interaction), collapse = ", ")) %>%
    filter(interactions == given_interaction) %>%
    pull(CHR_Group)
}

chr_group_each_interaction <- lapply(unique_interaction_vars, function(int) {
  find_chr_group_by_interaction(partition, int)
})

names(chr_group_each_interaction) <- unique_interaction_vars
all_chr_group = unlist(chr_group_each_interaction, use.names=FALSE)
print(length(all_chr_group))


# get swaooed data
swapped.genotypes.filename <-  paste0("swapped_genotypes/chr_merged/UKB_swp_", my_population, 
                                      "_res", resolution, "_ordered_whr")

cat("Attaching bigSNP object... ")
if(!file.exists(paste0(swapped.genotypes.filename, ".bk"))) {
  swapped.genotypes.filename = snp_readBed2(paste0(swapped.genotypes.filename, ".bed"))
}
obj.bigSNP.swapped <- snp_attach(paste0(swapped.genotypes.filename, ".rds"))
cat("done.\n")

G.swapped = obj.bigSNP.swapped$genotypes
map.swapped = obj.bigSNP.swapped$map %>%
  separate_wider_delim(marker.ID, ".", names = c("marker.ID", "knock_ind"))

# merge in group information
grp.file <- c()
for(c in seq(1, 22)) {
  myfilename <- paste0("group_information/ukb_gen_chr", c, "_ibd1_res", resolution, "_grp",".txt")
  myfile <- read.csv(myfilename, sep="")
  myfile$CHR <- c
  grp.file <- rbind(grp.file, myfile)
}

# merge group information into map
map.swapped = map.swapped %>%
  rename(SNP = marker.ID,
         CHR = chromosome) %>%
  left_join(grp.file, by = c("CHR", "SNP")) %>%
  mutate(CHR_Group = paste0(CHR, "_", Group))

map.swapped = map.swapped %>%
  mutate(Knockoff = ifelse(knock_ind == "a", FALSE, TRUE)) %>%
  mutate(CHR_Group_Knock = paste0(CHR_Group, "_", Knockoff))

# save information about loci within each partition
list_for_combine_lasso = list(unique_interaction_vars = unique_interaction_vars, 
                              all_environments = lapply(all_environments, function(x) sort(unique(x))))

write_rds(list_for_combine_lasso, file = paste0("lasso_each_env/combination_list_",
                                    my_population, "_res", resolution, "_", pheno.name, "_cov_", 
                                    covar_identifiyer, "_", suffix, "_", num.int))


# submit separate lasso jobs for each partition
for(i in 1:length(unique_interaction_vars)) {

  print(paste0("Working on interaction ", unique_interaction_vars[[i]]))
  interaction_character = paste(unique_interaction_vars[[i]], collapse = ".")
  
  # get all variables with this interaction combination
  relevant_chr_group = chr_group_each_interaction[[i]]


  # get environments for this interaction
  environments = all_environments[[i]]
  unique_environments = sort(unique(environments))


  for(e in unique_environments) {

    print(paste0("Working on environment ", e, " out of ", length(unique_environments)))

    index_e = which(environments == e)
    index_chr_group = which(map.swapped$CHR_Group %in% relevant_chr_group)
    
    index_list = list(environment_index_obs = index_e, 
                      relevant_chr_group_index = index_chr_group)
    
    write_rds(index_list, file = paste0("/lasso_each_env/",
                                          my_population, "_res", resolution, "_", pheno.name, "_cov_", covar_identifiyer, "_", 
                                        interaction_character, "_env", e, "_", suffix, "_", num.int))

    
    # submit shell file with lasso
    submit_lasso <- function() {

      # Create .sh file to submit jobs
      filename <- paste0("sh_files_env_lasso/submit_lasso_", 
                         my_population, "_",resolution, "_",pheno.name, "_cov_", covar_identifiyer,
                         "_", interaction_character, "_", e, "_", suffix , "_", num.int,".sh")

      # Open file connection for writing
      con <- file(filename, "w")

      # Write content to the file
      cat("#!/bin/bash\n", file = con)
      cat("#SBATCH --job-name=ls\n", file = con)
      cat(sprintf("#SBATCH --time=%d:00:00\n", time_lasso),file = con)
      cat(sprintf("#SBATCH --mem-per-cpu=%dGB\n", memory_lasso), file = con)
      cat(sprintf("#SBATCH --partition=%s\n", partition_sherlock), file = con)
      cat("#SBATCH --output=slurm_files/slurm-%j.out\n", file = con)
      cat("\n", file = con)
      cat("# Save job info on joblog:\n", file = con)
      cat("echo 'Job $JOB_ID started on:   ' `hostname -s`\n", file = con)
      cat("echo 'Job $JOB_ID started on:   ' `date `\n", file = con)
      cat("\n", file = con)
      cat("# Run code\n", file = con)
      cat("ml R/4.2.0\n", file = con)
      cat(sprintf("Rscript utils/run_lasso_function_environment.R %s %d %s %s %d %s %d %d %s %d %s\n",
                  my_population, resolution, pheno.name, 
                  interaction_character, e, pheno.file, dfmax, ncores, suffix, num.int, real_data_file), file = con)
      cat("\n", file = con)
      cat("# Echo job info on joblog:\n", file = con)
      cat("echo 'Job $JOB_ID ended on:   ' `hostname -s`\n", file = con)
      cat("echo 'Job $JOB_ID ended on:   ' `date `\n", file = con)
      cat("#echo ' '\n", file = con)


      # Close the file connection
      close(con)

      # Submit job using sbatch
      system(paste("sbatch", filename))

    }

    # Call the function
    submit_lasso()
    
    
  }

}

