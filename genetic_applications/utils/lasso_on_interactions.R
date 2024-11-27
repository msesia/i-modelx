#!/usr/bin/env Rscript

# UKB GWAS
#
# Class: R script
#
# This script runs the lasso on the genotypes with interactions


# Input arguments
args <- commandArgs(trailingOnly=TRUE)
my_population       <- as.character(args[1])
resolution       <- as.integer(args[2])
pheno.name       <- as.character(args[3])
dfmax       <- as.integer(args[4])
ncores       <- as.integer(args[5])
num.int  <- as.integer(args[6])
pheno.file <- as.character(args[7])

library(bigsnpr)
library(tidyverse)

set.seed(2024)

######## Define paths + filenames ################ 

covariate_file = paste0("swapped_genotypes/covariates/covariates_", pheno.name, "_",
                        my_population ,"_sample_all_groups")
covariates_interactions = read.csv(covariate_file)

# get last part after "_" for covariate naming
get_last_string <- function(string) {
  parts <- strsplit(string, "_")[[1]]
  return(tail(parts, 1))
}
filtered_col_names <- colnames(covariates_interactions)[!grepl("eid", colnames(covariates_interactions)) & !grepl("^rev", colnames(covariates_interactions))]
covar_names_pretty <- sapply(filtered_col_names, get_last_string)
covar_identifiyer = paste0(sort(covar_names_pretty), collapse = "_")


obj.bigSNP.swapped.interaction.filename = paste0("swapped_genotypes/covariate_interactions/merged/UKB_swp_",
                                                 my_population, "_res", resolution, "_", pheno.name,"_covar_interaction_", covar_identifiyer)

stats.path = paste0("UKB_stats/", 
                    my_population, "_res", resolution, "_", pheno.name, "_cov_", covar_identifiyer, "_numint", num.int)
                    

######### Load file with all interactions ###########


if(!file.exists(paste0(obj.bigSNP.swapped.interaction.filename, ".bk"))) {
  obj.bigSNP.swapped.interaction.rds = snp_readBed2(paste0(obj.bigSNP.swapped.interaction.filename, ".bed"))
}
cat("done.\n")

obj.bigSNP.swapped.interaction <- snp_attach(paste0(obj.bigSNP.swapped.interaction.filename, ".rds"))

G.swapped.interaction = obj.bigSNP.swapped.interaction$genotypes

# Extract list of variants
map <- obj.bigSNP.swapped.interaction$map %>% as_tibble()
colnames(map) <- c("CHR", "SNP", "gd", "BP", "a1", "a2")
map <- map %>% 
  select(CHR, SNP, BP) %>% 
  separate_wider_delim(SNP, delim = ".", names = c("SNP", "swapped_knockoff")) %>% 
  separate_wider_delim(swapped_knockoff, delim = ":", names = c("swapped_knockoff", "Interaction"), too_few = "align_start") %>% 
  mutate(Interaction = ifelse(is.na(Interaction), "none", Interaction), 
         Knockoff = ifelse(swapped_knockoff == "b", TRUE, FALSE))

# merge in group information
grp.file <- c()
for(c in seq(1, 22)) {
  myfilename <- paste0("group_information/ukb_gen_chr", c, "_ibd1_res", resolution, "_grp",".txt")
  myfile <- read.csv(myfilename, sep="") 
  myfile$CHR <- c
  grp.file <- rbind(grp.file, myfile)
}

map = map %>% left_join(grp.file, by = c("CHR", "SNP"))

# Extract list of subjects
Subjects <- obj.bigSNP.swapped.interaction$fam %>% as_tibble()
colnames(Subjects) <- c("FID", "IID", "X1", "X2", "sex", "X3")
Subjects <- Subjects %>% select(FID, IID) %>%
  mutate(FID=as.character(FID), IID=as.character(IID))

############# Scaling factors, covariates, phenotypes ###############

# Compute scaling factor for the genotypes
cat("Computing scaling factors for all ... ")
scaler <- big_scale()
G.scale <- scaler(G.swapped.interaction)
scaling.factors <- G.scale$scale
cat("done.\n")

cat("Reading phenotype file... ")
Phenotypes <- read.csv(pheno.file)
cat("done.\n")

if(sum(Phenotypes$eid != Subjects$FID) > 0) {
  stop("Subjects not in correct order!")
}

# Extract response variable
ind.train <- which(!is.na(Phenotypes[[pheno.name]]))
cat(sprintf("%d individuals out of %d have missing phenotype.\n",
            nrow(Subjects)-length(ind.train), nrow(Subjects)))
y <- Phenotypes[[pheno.name]][ind.train]


if(sum(covariates_interactions$eid != Subjects$FID) > 0) {
  stop("Subjects not in correct order!")
}

# get covariates
baseline_covariate_names = c("age", "age_squared", "sex","PC1", "PC2", "PC3", "PC4", "PC5", "eid")

names_covariates_interactions = colnames(covariates_interactions %>% dplyr::select(-eid))

unique_baseline_covariate_names = baseline_covariate_names[!(baseline_covariate_names %in% names_covariates_interactions)] 
baseline_covariates <- Phenotypes %>% dplyr::select(all_of(unique_baseline_covariate_names))

covar.train = baseline_covariates %>% left_join(covariates_interactions, by = "eid") %>% dplyr::select(-eid)
covar.train = covar.train %>% dplyr::select(-starts_with("rev"))
covar.train = as.matrix(covar.train)[ind.train,,drop=F] 

########### Running the lasso + get coefficients ########

weight_penalty_vector = rep(1, ncol(G.swapped.interaction))
weight_penalty_vector[which(map$Interaction == "none")] <- 0.25

cat(sprintf("Fitting sparse linear regression with %d observations, %d variants and %d covariates... ",
            length(y), ncol(G.swapped.interaction), ncol(covar.train)))
lasso.fit <- big_spLinReg(G.swapped.interaction, y.train=y, ind.train=ind.train, covar.train=covar.train,
                          dfmax=dfmax, ncores=ncores, pf.X = weight_penalty_vector) 
cat("done.\n")


# Extract beta from each fold and combine them
cat("Extracting regression coefficients... ")
beta <- sapply(1:10, function(k) lasso.fit[[1]][k][[1]]$beta)

# Separate the coefficients of the genetic variants from the coefficients of the covariates
total_number_coefficients = nrow(beta)
total_number_genotypes_in_model = total_number_coefficients - ncol(covar.train)
beta.variants <- beta[1:total_number_genotypes_in_model,]
beta.covariates <- beta[(total_number_genotypes_in_model+1):nrow(beta),]

# Undo scaling of lasso coefficients
remaining_columns = attr(lasso.fit, "ind.col") # which genotypes were used for model fitting
beta.variants <- beta.variants * scaling.factors[remaining_columns]
Beta <- cbind(tibble(CHR=map[remaining_columns, ]$CHR,
                     SNP=map[remaining_columns, ]$SNP,
                     BP=map[remaining_columns, ]$BP, 
                     Knockoff=map[remaining_columns, ]$Knockoff, 
                     Interaction = map[remaining_columns, ]$Interaction),
              as_tibble(beta.variants)) %>% as_tibble()
colnames(Beta) <- c("CHR", "SNP", "BP", "Knockoff", "Interaction",paste("K", seq(ncol(beta.variants)),sep=""))
Beta <- Beta %>%
  mutate(Z=(K1+K2+K3+K4+K5+K6+K7+K8+K9+K10)/10,
         Nonzero=(K1!=0)+(K2!=0)+(K3!=0)+(K4!=0)+(K5!=0)+(K6!=0)+(K7!=0)+(K8!=0)+(K9!=0)+(K10!=0)) %>%
  select(CHR, SNP, BP, Knockoff, Interaction,Z)

# Extract the estimated coefficients
Lasso.res <- Beta %>%
  inner_join(map, by = c("CHR", "SNP", "BP", "Knockoff", "Interaction")) %>%
  select(CHR, SNP, BP, Z, Group, Knockoff, Interaction)
cat("done.\n")

################## get top interactions ################

W.stats.sum <- function(Z, knockoff) {
  importance <- abs(Z)
  z  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
  zk <- sum(importance[which(knockoff==TRUE)], na.rm=T)
  w <- z+zk
}

Stats <- Lasso.res %>%
  select("CHR", "Group", "SNP", "BP", "Knockoff", "Z", "Interaction") %>%
  group_by(CHR, Group, Interaction) %>%
  summarize(W.sum = W.stats.sum(abs(Z),Knockoff),
            Lead=which.max(abs(Z)), SNP.lead=SNP[Lead], BP.lead=BP[Lead],
            Size=n()) %>%
  ungroup() %>%
  select(CHR, Group, SNP.lead, BP.lead, Size, W.sum, Interaction)
cat("done.\n")

Stats = Stats %>% 
  rename(Interaction_orig = Interaction) %>%
  mutate(Interaction = sub("^rev_", "", Interaction_orig))

partition_tmp <- Stats %>%
  filter(W.sum != 0) %>%
  group_by(CHR, Group) %>% 
  dplyr::mutate(number_unique_interactions = n_distinct(Interaction))  %>% 
  arrange(CHR, Group) %>% 
  dplyr::filter(!(number_unique_interactions > 1 & Interaction == "none")) %>% 
  # if there are two interactions picked for the same covariate (e.g. with sex and rev_sex) pick only the top one
  group_by(CHR, Group, Interaction) %>% 
  top_n(pmin(n(),1), W.sum) 
  
partition = partition_tmp %>% 
  group_by(CHR, Group) %>%
  top_n(pmin(n(),num.int), W.sum) %>% 
  dplyr::select(-number_unique_interactions)


# Keep an interaction only if it has a main effect as asll
chr_group_with_main_effect = Stats %>% filter(W.sum != 0) %>%
  filter(Interaction == "none") %>% 
  dplyr::select(CHR, Group) %>% 
  unique() %>% 
  mutate(has_main_effect = 1)

partition_with_main_tmp <- Stats %>%
  left_join(chr_group_with_main_effect, by = c("CHR", "Group")) %>% 
  filter(has_main_effect == 1) %>%
  filter(W.sum != 0) %>%
  group_by(CHR, Group) %>% 
  dplyr::mutate(number_unique_interactions = n_distinct(Interaction))  %>% 
  arrange(CHR, Group) %>% 
  dplyr::filter(!(number_unique_interactions > 1 & Interaction == "none")) %>% 
  # if there are two interactions picked for the same covariate (e.g. with sex and rev_sex) pick only the top one
  group_by(CHR, Group, Interaction) %>% 
  top_n(pmin(n(),1), W.sum) 

partition_with_main = partition_with_main_tmp %>%
  group_by(CHR, Group) %>%
  top_n(pmin(n(),num.int), W.sum) %>% 
  dplyr::select(-number_unique_interactions)

print("rows partition:")
print(nrow(partition))

print("rows partition (restricting to only those with main effect)")
print(nrow(partition_with_main))

print("Average W.sum for main effects")
round(summary(Stats %>% filter(W.sum != 0) %>% filter(Interaction == "none") %>% pull(W.sum)), 5)

print("Average W.sum in partition for interaction effects (if also having main effect)")
round(summary(partition_with_main %>% ungroup %>% filter(Interaction != "none") %>% pull(W.sum)), 5)

print("Average W.sum in partition for interaction effects (no requirement on main effect)")
round(summary(partition %>% ungroup %>% filter(Interaction != "none") %>% pull(W.sum)), 5)

# alternative partition: interaction only if interaction W is greater than 2 * median effect
median_main_effect = median(Stats %>% filter(W.sum != 0) %>% filter(Interaction == "none") %>% pull(W.sum))

partition_interaction_only_if_larger_med_main_tmp <- Stats %>%
  filter(W.sum != 0) %>%
  filter(!(Interaction != "none" & W.sum < 2*median_main_effect)) %>%
  mutate(less_than_median = (W.sum < median_main_effect) & Interaction != "none") %>% 
  group_by(CHR, Group) %>% #maybe change later: choose by chr, group, SNP.raw
  dplyr::mutate(number_unique_interactions = n_distinct(Interaction))  %>% 
  arrange(CHR, Group) %>% 
  dplyr::filter(!(number_unique_interactions > 1 & Interaction == "none")) %>% 
  # if there are two interactions picked (e.g. with sex and rev_sex) pick only the top one
  group_by(CHR, Group, Interaction) %>% 
  top_n(pmin(n(),1), W.sum) 

partition_interaction_only_if_larger_med_main = partition_interaction_only_if_larger_med_main_tmp %>% 
  group_by(CHR, Group) %>%
  top_n(pmin(n(),num.int), W.sum) %>% 
  dplyr::select(-number_unique_interactions)

print("Average W.sum in partition")
round(summary(partition_interaction_only_if_larger_med_main %>% ungroup %>% filter(Interaction != "none") %>% pull(W.sum)), 5)

print("Average W.sum in main effect")
round(summary(partition_interaction_only_if_larger_med_main %>% ungroup %>% filter(Interaction == "none") %>% pull(W.sum)), 5)

# alternative partition: interaction only if interaction effect is at least twice as large as corresponding main effect
main_effect_sizes = Stats %>% 
  filter(Interaction == "none") %>% 
  select(CHR, Group, W.sum) %>% 
  rename(W.sum.main = W.sum)

partition_interaction_only_if_larger_main_tmp <- Stats %>%
  filter(W.sum != 0) %>%
  left_join(main_effect_sizes, by = c("CHR", "Group")) %>%
  filter(!(Interaction != "none" & W.sum < 2*W.sum.main)) %>%
  mutate(less_than_median = (W.sum < median_main_effect) & Interaction != "none") %>% 
  group_by(CHR, Group) %>% #maybe change later: choose by chr, group, SNP.raw
  dplyr::mutate(number_unique_interactions = n_distinct(Interaction))  %>% 
  arrange(CHR, Group) %>% 
  dplyr::filter(!(number_unique_interactions > 1 & Interaction == "none")) %>% 
  # if there are two interactions picked (e.g. with sex and rev_sex) pick only the top one
  group_by(CHR, Group, Interaction) %>% 
  top_n(pmin(n(),1), W.sum) 

partition_interaction_only_if_larger_main = partition_interaction_only_if_larger_main_tmp %>% 
  group_by(CHR, Group) %>%
  top_n(pmin(n(),num.int), W.sum) %>% 
  dplyr::select(-number_unique_interactions)

# alternative partition: keep main effect for all, also for those that are 0
partition_all_main_tmp <- Stats %>%
  filter(!(W.sum == 0 & Interaction != "none")) %>% # keep main effects, regardless of whether W = 0 or not, only drop zero interactions
  filter(!(Interaction != "none" & W.sum < median_main_effect)) %>%
  group_by(CHR, Group) %>% 
  dplyr::mutate(number_unique_interactions = n_distinct(Interaction))  %>% 
  arrange(CHR, Group) %>% 
  dplyr::filter(!(number_unique_interactions > 1 & Interaction == "none")) %>% 
  # if there are two interactions picked (e.g. with sex and rev_sex) pick only the top one
  group_by(CHR, Group, Interaction) %>% 
  top_n(pmin(n(),1), W.sum) 

partition_all_main = partition_all_main_tmp %>% 
  group_by(CHR, Group) %>%
  top_n(pmin(n(),num.int), W.sum) %>% 
  dplyr::select(-number_unique_interactions)

############ Save results ########## 

# save SNPs in Groups that are remaining (that are nonzero)
# only unique SNPs in the group
nonzero_snps_in_group = Lasso.res %>% 
  dplyr::select(CHR, Group, SNP, BP) %>% 
  unique() %>% 
  dplyr::filter(CHR %in% partition$CHR & Group %in% partition$Group)

# Save results
covariate_initials = paste0(gsub("[^a-zA-Z]", "", gsub("([^_])[^_]*", "\\1", colnames(covariates_interactions))), collapse = "")

stats.file <- sprintf("%s_lasso_interaction_stats.txt", stats.path)
Stats %>% write_delim(stats.file, delim=" ")
cat(sprintf("Test statistics written on:\n%s\n", stats.file))

partition.file <- sprintf("%s_lasso_interaction_partition.txt", stats.path)
partition %>% write_delim(partition.file, delim=" ")
cat(sprintf("Partition written on:\n%s\n", partition.file))

partition.file.main <- sprintf("%s_lasso_interaction_partition_main.txt", stats.path)
partition_with_main %>% write_delim(partition.file.main, delim=" ")
cat(sprintf("Partition main effects written on:\n%s\n", partition.file.main))

partition.file.filtered <- sprintf("%s_lasso_interaction_partition_filtered.txt", stats.path)
partition_interaction_only_if_larger_med_main %>% write_delim(partition.file.filtered, delim=" ")
cat(sprintf("Partition main effects written on:\n%s\n", partition.file.filtered))

partition.file.filtered.larger <- sprintf("%s_lasso_interaction_partition_filtered_larger.txt", stats.path)
partition_interaction_only_if_larger_main %>% write_delim(partition.file.filtered.larger, delim=" ")
cat(sprintf("Partition main effects written on:\n%s\n", partition.file.filtered.larger))

partition.file.all.main <- sprintf("%s_lasso_interaction_partition_all_main.txt", stats.path)
partition_all_main %>% write_delim(partition.file.all.main, delim=" ")
cat(sprintf("Partition main effects written on:\n%s\n", partition.file.all.main))

nonzero.snp.file <- sprintf("%s_lasso_interaction_nonzero_snps.txt", stats.path)
nonzero_snps_in_group %>% write_delim(nonzero.snp.file, delim=" ")
cat(sprintf("Nonzero snps written on:\n%s\n", nonzero.snp.file))

