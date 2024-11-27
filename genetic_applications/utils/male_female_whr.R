#!/usr/bin/env Rscript

# UK Biobank GWAS
#
# Class: R script
#
# This script runs the knockoff filter separately for males and females (for WHR) as well as for the entire sample


# Input arguments
args <- commandArgs(trailingOnly=TRUE)
my_population       <- as.character(args[1])
resolution       <- as.integer(args[2])
pheno.name       <- as.character(args[3])
fdr  <- as.numeric(args[4])
dfmax <- as.numeric(args[5])
ncores <- as.numeric(args[6])
pheno.file       <- as.character(args[7])
real_data_file       <- as.character(args[8])


library(bigsnpr)
library(tidyverse)

set.seed(2022)

################## Define paths ################ 


# get covariates
covariate_file = paste0("covariates/covariates_", pheno.name, "_",
                        my_population ,"_sample_all_groups")

print(covariate_file)

covariates_interactions = read.csv(covariate_file)
head(covariates_interactions)
names_covariates_interactions = colnames(covariates_interactions %>% dplyr::select(-eid))

print(summary(covariates_interactions))

get_last_string <- function(string) {
  parts <- strsplit(string, "_")[[1]]
  return(tail(parts, 1))
}
filtered_col_names <- colnames(covariates_interactions)[!grepl("eid", colnames(covariates_interactions)) & !grepl("^rev", colnames(covariates_interactions))]
covar_names_pretty <- sapply(filtered_col_names, get_last_string)
covar_identifiyer = paste0(sort(covar_names_pretty), collapse = "_")

print("Covar identifier")
print(covar_identifiyer)

# path where to save stats
stats.path = paste0("UKB_stats/", my_population, 
                    "_res", resolution, "_", pheno.name, "_cov_", covar_identifiyer)

W.stats <- function(Z, knockoff) {
  importance <- abs(Z)
  z  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
  zk <- sum(importance[which(knockoff==TRUE)], na.rm=T)
  w <- z-zk
}

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


################ get real data ############

bed_bim_fam_filename <- real_data_file

# Attach the "bigSNP" object in R session
cat("Attaching bigSNP object... ")
obj.bigSNP.real.full <- snp_attach(paste0(bed_bim_fam_filename,".rds"))
cat("done.\n")

######### get phenotype ###########

cat("Reading phenotype file... ")
Phenotypes <- read.csv(pheno.file)
cat("done.\n")

ind.train <- which(!is.na(Phenotypes[[pheno.name]]))
y <- Phenotypes[[pheno.name]][ind.train]


# include all covariates, including those that we are looking for interactions for
baseline_covariate_names = c("age", "age_squared", "sex","PC1", "PC2", "PC3", "PC4", "PC5", "eid")
names_covariates_interactions = colnames(covariates_interactions %>% dplyr::select(-eid))

print(names_covariates_interactions)
print(baseline_covariate_names)

unique_baseline_covariate_names = baseline_covariate_names[!(baseline_covariate_names %in% names_covariates_interactions)] 
print(unique_baseline_covariate_names)

baseline_covariates <- Phenotypes %>% dplyr::select(all_of(unique_baseline_covariate_names))

head(baseline_covariates)
head(covariates_interactions)
covar.train = baseline_covariates %>% dplyr::left_join(covariates_interactions, by = "eid") %>% dplyr::select(-eid)
covar.train = covar.train %>% dplyr::select(-starts_with("rev"))
covar.train = as.matrix(covar.train)[ind.train,,drop=F] 

print(summary(covar.train))

################## for the entire sample ###################

G.real = obj.bigSNP.real.full$genotypes

print("dimensions entire genotype")
print(dim(G.real))

# get list of variants
map.real = obj.bigSNP.real.full$map %>% as_tibble()
colnames(map.real) <- c("CHR", "SNP", "gd", "BP", "a1", "a2")
map.real <- map.real %>% dplyr::select(CHR, SNP, BP) %>%
  mutate(SNP_knock = SNP,
         Knockoff = ifelse(endsWith(SNP, ".k"), TRUE, FALSE),
         SNP = sub(".k$", "", SNP)) %>%
  dplyr::select(CHR, SNP, SNP_knock, BP, Knockoff)

# merge in group information
grp.file <- c()
chr <- seq(1, 22)
for(c in chr) {
  myfilename <- paste0("group_information/ukb_gen_chr", c, "_ibd1_res", resolution, "_grp",".txt")
  myfile <- read.csv(myfilename, sep="")
  
  myfile$CHR <- c
  
  grp.file <- rbind(grp.file, myfile)
}

map.real = map.real %>% left_join(grp.file, by = c("CHR", "SNP")) %>%
  mutate(CHR_Group = paste0(CHR, "_", Group))  %>%
  mutate(CHR_Group_Knock = paste0(CHR_Group, "_", Knockoff))

# run lasso 
cat("Run lasso directly on real data ... \n")
scaler <- big_scale()
G.scale <- scaler(G.real)
scaling.factors <- G.scale$scale

# Run lasso
cat("Running lasso... ")
lasso.fit <- big_spLinReg(G.real,
                          y.train=y,
                          dfmax=dfmax,
                          ncores=ncores,
                          covar.train = covar.train,
                          pf.covar=rep(0, ncol(covar.train)))

# Extract beta from each fold and combine them
beta <- sapply(1:10, function(k) lasso.fit[[1]][k][[1]]$beta)

# Separate the coefficients of the genetic variants from the coefficients of the covariates
total_number_coefficients = nrow(beta)
total_number_genotypes_in_model = total_number_coefficients - ncol(covar.train)
beta.variants <- beta[1:total_number_genotypes_in_model,]

# Undo scaling of lasso coefficients
remaining_columns = attr(lasso.fit, "ind.col") # which genotypes were used for model fitting (index does not include covariates)
beta.variants <- beta.variants * scaling.factors[remaining_columns]
Beta <- cbind(tibble(CHR=map.real[remaining_columns, ]$CHR,
                     SNP=map.real[remaining_columns, ]$SNP,
                     BP=map.real[remaining_columns, ]$BP,
                     Knockoff=map.real[remaining_columns, ]$Knockoff),
              as_tibble(beta.variants)) %>% as_tibble()
colnames(Beta) <- c("CHR", "SNP", "BP", "Knockoff",paste("K", seq(ncol(beta.variants)),sep=""))
Beta <- Beta %>%
  mutate(Z=(K1+K2+K3+K4+K5+K6+K7+K8+K9+K10)/10,
         Nonzero=(K1!=0)+(K2!=0)+(K3!=0)+(K4!=0)+(K5!=0)+(K6!=0)+(K7!=0)+(K8!=0)+(K9!=0)+(K10!=0)) %>%
  select(CHR, SNP, BP, Knockoff,Z)

Stats_real <- Beta %>%
  inner_join(map.real, by = c("CHR", "SNP", "BP", "Knockoff")) %>%
  select("CHR", "Group", "SNP", "BP", "Knockoff", "Z") %>%
  group_by(CHR, Group) %>%
  summarize(W = W.stats(abs(Z),Knockoff),
            Lead=which.max(abs(Z)), SNP.lead=SNP[Lead], BP.lead=BP[Lead],
            Size=n()) %>%
  ungroup() %>%
  select(CHR, Group, SNP.lead, BP.lead, Size, W)


selections_real <- Stats_real %>% knockoff.filter(fdr=fdr, offset=1) %>%
  arrange(CHR, Group) %>%
  arrange(desc(W))

print(paste0("Number rejections directly on real data:", nrow(selections_real)))


# # Save list of discoveries
selections.out.file <- sprintf("%s_lasso_discoveries_real_mf.txt", stats.path)
selections_real %>% write_delim(selections.out.file, delim=" ")
cat(sprintf("lasso discoveries female real directly: %s\n", selections.out.file))

stats.out.file <- sprintf("%s_lasso_stats_real_mf.txt", stats.path)
Stats_real %>% write_delim(stats.out.file, delim=" ")
cat(sprintf("lasso stats female real directly written on: %s\n", stats.out.file))


################## SPLIT INTO MALE AND FEMALE ########################## 

female_not_missing_indices <- intersect(ind.train, which(Phenotypes$sex == 0))
male_not_missing_indices <- intersect(ind.train, which(Phenotypes$sex == 1))

# save male and female data
backingfile_name_female = paste0("female_male/ukb_",
                                 my_population,
                                 "_female",
                                 "_res",
                                 resolution)

if(file.exists(paste0(backingfile_name_female, ".bk"))){
  file.remove(paste0(backingfile_name_female, ".bk"))
  file.remove(paste0(backingfile_name_female, ".rds"))
  file.remove(paste0(backingfile_name_female, ".bed"))
  file.remove(paste0(backingfile_name_female, ".bim"))
  file.remove(paste0(backingfile_name_female, ".fam"))
}


obj.bigSNP.real.female = snp_subset(obj.bigSNP.real.full, 
                                    ind.row = female_not_missing_indices,
                                    backingfile = backingfile_name_female)

backingfile_name_male = paste0("female_male/ukb_",
                               my_population,
                               "_male",
                               "_res",
                               resolution)

if(file.exists(paste0(backingfile_name_male, ".bk"))){
  file.remove(paste0(backingfile_name_male, ".bk"))
  file.remove(paste0(backingfile_name_male, ".rds"))
  file.remove(paste0(backingfile_name_male, ".bed"))
  file.remove(paste0(backingfile_name_male, ".bim"))
  file.remove(paste0(backingfile_name_male, ".fam"))
}


obj.bigSNP.real.male = snp_subset(obj.bigSNP.real.full, 
                                  ind.row = male_not_missing_indices,
                                  backingfile = backingfile_name_male)

##### ___For females #######


cat("Attaching bigSNP object for females ... ")
obj.bigSNP.real.female <- snp_attach(paste0(obj.bigSNP.real.female))
cat("done.\n")


y.female <- Phenotypes[[pheno.name]][female_not_missing_indices]

covar.train.female = baseline_covariates %>% left_join(covariates_interactions, by = "eid") %>% dplyr::select(-eid)
covar.train.female = covar.train.female %>% dplyr::select(-starts_with("rev")) %>% dplyr::select(-sex)
covar.train.female = as.matrix(covar.train.female)[female_not_missing_indices,,drop=F] 

G.real.female = obj.bigSNP.real.female$genotypes

print("dimensions female genotype")
print(dim(G.real.female))

# get list of variants
map.real.female = obj.bigSNP.real.female$map %>% as_tibble()
colnames(map.real.female) <- c("CHR", "SNP", "gd", "BP", "a1", "a2")
map.real.female <- map.real.female %>% dplyr::select(CHR, SNP, BP) %>%
  mutate(SNP_knock = SNP,
         Knockoff = ifelse(endsWith(SNP, ".k"), TRUE, FALSE),
         SNP = sub(".k$", "", SNP)) %>%
  dplyr::select(CHR, SNP, SNP_knock, BP, Knockoff)

map.real.female = map.real.female %>% left_join(grp.file, by = c("CHR", "SNP")) %>%
  mutate(CHR_Group = paste0(CHR, "_", Group))  %>%
  mutate(CHR_Group_Knock = paste0(CHR_Group, "_", Knockoff))

# run lasso 
cat("Run lasso directly on real data ... \n")
scaler <- big_scale()
G.scale.female <- scaler(G.real.female)
scaling.factors.female <- G.scale.female$scale

# Run lasso
cat("Running lasso... ")
lasso.fit.female <- big_spLinReg(G.real.female,
                                 y.train=y.female,
                                 dfmax=dfmax,
                                 ncores=ncores,
                                 covar.train = covar.train.female,
                                 pf.covar=rep(0, ncol(covar.train.female)))

# Extract beta from each fold and combine them
beta.female <- sapply(1:10, function(k) lasso.fit.female[[1]][k][[1]]$beta)

# Separate the coefficients of the genetic variants from the coefficients of the covariates
total_number_coefficients_female = nrow(beta.female)
total_number_genotypes_in_model_female = total_number_coefficients_female - ncol(covar.train.female)
beta.variants.female <- beta.female[1:total_number_genotypes_in_model_female,]

# Undo scaling of lasso coefficients
remaining_columns_female = attr(lasso.fit.female, "ind.col") # which genotypes were used for model fitting (index does not include covariates)
beta.variants.female <- beta.variants.female * scaling.factors.female[remaining_columns_female]
Beta.female <- cbind(tibble(CHR=map.real.female[remaining_columns_female, ]$CHR,
                            SNP=map.real.female[remaining_columns_female, ]$SNP,
                            BP=map.real.female[remaining_columns_female, ]$BP,
                            Knockoff=map.real.female[remaining_columns_female, ]$Knockoff),
                     as_tibble(beta.variants.female)) %>% as_tibble()
colnames(Beta.female) <- c("CHR", "SNP", "BP", "Knockoff",paste("K", seq(ncol(beta.variants.female)),sep=""))
Beta.female <- Beta.female %>%
  mutate(Z=(K1+K2+K3+K4+K5+K6+K7+K8+K9+K10)/10,
         Nonzero=(K1!=0)+(K2!=0)+(K3!=0)+(K4!=0)+(K5!=0)+(K6!=0)+(K7!=0)+(K8!=0)+(K9!=0)+(K10!=0)) %>%
  select(CHR, SNP, BP, Knockoff,Z)

Stats_real_female <- Beta.female %>%
  inner_join(map.real.female, by = c("CHR", "SNP", "BP", "Knockoff")) %>%
  select("CHR", "Group", "SNP", "BP", "Knockoff", "Z") %>%
  group_by(CHR, Group) %>%
  summarize(W = W.stats(abs(Z),Knockoff),
            Lead=which.max(abs(Z)), SNP.lead=SNP[Lead], BP.lead=BP[Lead],
            Size=n()) %>%
  ungroup() %>%
  select(CHR, Group, SNP.lead, BP.lead, Size, W)



selections_real_female <- Stats_real_female %>% knockoff.filter(fdr=fdr, offset=1) %>%
  arrange(CHR, Group) %>%
  arrange(desc(W))

print(paste0("Number rejections female directly on real data:", nrow(selections_real_female)))


# Save list of discoveries
selections.out.file <- sprintf("%s_lasso_discoveries_female_real.txt", stats.path)
selections_real_female %>% write_delim(selections.out.file, delim=" ")
cat(sprintf("lasso discoveries female real directly: %s\n", selections.out.file))

stats.out.file <- sprintf("%s_lasso_stats_female_real.txt", stats.path)
Stats_real_female %>% write_delim(stats.out.file, delim=" ")
cat(sprintf("lasso stats female real directly written on: %s\n", stats.out.file))



##### ___For males #######

cat("Attaching bigSNP object for males ... ")
obj.bigSNP.real.male <- snp_attach(paste0(obj.bigSNP.real.male))
cat("done.\n")

y.male <- Phenotypes[[pheno.name]][male_not_missing_indices]

covar.train.male = baseline_covariates %>% left_join(covariates_interactions, by = "eid") %>% dplyr::select(-eid)
covar.train.male = covar.train.male %>% dplyr::select(-starts_with("rev")) %>% dplyr::select(-sex)
covar.train.male = as.matrix(covar.train.male)[male_not_missing_indices,,drop=F] 

G.real.male = obj.bigSNP.real.male$genotypes

print("dimensions male genotype")
print(dim(G.real.male))

# get list of variants
map.real.male = obj.bigSNP.real.male$map %>% as_tibble()
colnames(map.real.male) <- c("CHR", "SNP", "gd", "BP", "a1", "a2")
map.real.male <- map.real.male %>% dplyr::select(CHR, SNP, BP) %>%
  mutate(SNP_knock = SNP,
         Knockoff = ifelse(endsWith(SNP, ".k"), TRUE, FALSE),
         SNP = sub(".k$", "", SNP)) %>%
  dplyr::select(CHR, SNP, SNP_knock, BP, Knockoff)

map.real.male = map.real.male %>% left_join(grp.file, by = c("CHR", "SNP")) %>%
  mutate(CHR_Group = paste0(CHR, "_", Group))  %>%
  mutate(CHR_Group_Knock = paste0(CHR_Group, "_", Knockoff))

# run lasso 
cat("Run lasso directly on real data ... \n")
scaler <- big_scale()
G.scale.male <- scaler(G.real.male)
scaling.factors.male <- G.scale.male$scale

# Run lasso
cat("Running lasso... ")
lasso.fit.male <- big_spLinReg(G.real.male,
                               y.train=y.male,
                               dfmax=dfmax,
                               ncores=ncores,
                               covar.train = covar.train.male,
                               pf.covar=rep(0, ncol(covar.train.male)))

# Extract beta from each fold and combine them
beta.male <- sapply(1:10, function(k) lasso.fit.male[[1]][k][[1]]$beta)

# Separate the coefficients of the genetic variants from the coefficients of the covariates
total_number_coefficients_male = nrow(beta.male)
total_number_genotypes_in_model_male = total_number_coefficients_male - ncol(covar.train.male)
beta.variants.male <- beta.male[1:total_number_genotypes_in_model_male,]

# Undo scaling of lasso coefficients
remaining_columns_male = attr(lasso.fit.male, "ind.col") # which genotypes were used for model fitting (index does not include covariates)
beta.variants.male <- beta.variants.male * scaling.factors.male[remaining_columns_male]
Beta.male <- cbind(tibble(CHR=map.real.male[remaining_columns_male, ]$CHR,
                          SNP=map.real.male[remaining_columns_male, ]$SNP,
                          BP=map.real.male[remaining_columns_male, ]$BP,
                          Knockoff=map.real.male[remaining_columns_male, ]$Knockoff),
                   as_tibble(beta.variants.male)) %>% as_tibble()
colnames(Beta.male) <- c("CHR", "SNP", "BP", "Knockoff",paste("K", seq(ncol(beta.variants.male)),sep=""))
Beta.male <- Beta.male %>%
  mutate(Z=(K1+K2+K3+K4+K5+K6+K7+K8+K9+K10)/10,
         Nonzero=(K1!=0)+(K2!=0)+(K3!=0)+(K4!=0)+(K5!=0)+(K6!=0)+(K7!=0)+(K8!=0)+(K9!=0)+(K10!=0)) %>%
  select(CHR, SNP, BP, Knockoff,Z)

Stats_real_male <- Beta.male %>%
  inner_join(map.real.male, by = c("CHR", "SNP", "BP", "Knockoff")) %>%
  select("CHR", "Group", "SNP", "BP", "Knockoff", "Z") %>%
  group_by(CHR, Group) %>%
  summarize(W = W.stats(abs(Z),Knockoff),
            Lead=which.max(abs(Z)), SNP.lead=SNP[Lead], BP.lead=BP[Lead],
            Size=n()) %>%
  ungroup() %>%
  select(CHR, Group, SNP.lead, BP.lead, Size, W)

selections_real_male <- Stats_real_male %>% knockoff.filter(fdr=fdr, offset=1) %>%
  arrange(CHR, Group) %>%
  arrange(desc(W))

print(paste0("Number rejections male directly on real data:", nrow(selections_real_male)))

# Save list of discoveries
selections.out.file <- sprintf("%s_lasso_discoveries_male_real.txt", stats.path)
selections_real_male %>% write_delim(selections.out.file, delim=" ")
cat(sprintf("lasso discoveries male real directly: %s\n", selections.out.file))

stats.out.file <- sprintf("%s_lasso_stats_male_real.txt", stats.path)
Stats_real_male %>% write_delim(stats.out.file, delim=" ")
cat(sprintf("lasso stats male real directly written on: %s\n", stats.out.file))



######### ___Combined males + females ############

# concatenate the statistics from males and females and jointly pass them through the 
# knockoff filter

Stats_real_male = Stats_real_male %>% mutate(sample = "male") 
Stats_real_female = Stats_real_female %>% mutate(sample = "female") 

all_stats = rbind(Stats_real_male, Stats_real_female) %>%
  mutate(CHR_Group = paste0(CHR, "_", Group))

selections_jointly <- all_stats %>% knockoff.filter(fdr=fdr, offset=1) %>%
  arrange(CHR, Group) %>%
  arrange(desc(W))

print(paste0("Number rejections male/female directly on real data:", nrow(selections_jointly)))


# Save list of discoveries
selections.out.file <- sprintf("%s_lasso_discoveries_male_female_joint_real.txt", stats.path)
selections_jointly %>% write_delim(selections.out.file, delim=" ")
cat(sprintf("lasso discoveries male real directly: %s\n", selections.out.file))

stats.out.file <- sprintf("%s_lasso_stats_male_female_joint_real.txt", stats.path)
all_stats %>% write_delim(stats.out.file, delim=" ")
cat(sprintf("lasso stats male real directly written on: %s\n", stats.out.file))

# rejections that are made only a single time in either males or females 
joint_rejected_single_time = selections_jointly %>% 
  group_by(CHR_Group) %>% 
  mutate(num_rejected = n()) %>% 
  filter(num_rejected == 1)

print(paste0("Number rejected in joint analysis a single time: ", nrow(joint_rejected_single_time)))
print(paste0("Out of these: rejected female only: ", sum(joint_rejected_single_time$sample == "female")))
print(paste0("Out of these: rejected male only: ", sum(joint_rejected_single_time$sample == "male")))

