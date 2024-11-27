#!/usr/bin/env Rscript


# UKB GWAS
#
# Class: R script
#
# This script runs the vanilla knockoff on the original data and the pre-screened data

# Input arguments
args <- commandArgs(trailingOnly=TRUE)
my_population       <- as.character(args[1])
resolution       <- as.integer(args[2])
pheno.name       <- as.character(args[3])
fdr  <- as.numeric(args[4])
dfmax <- as.numeric(args[5])
ncores <- as.numeric(args[6])
pheno.file       <- as.character(args[7])
real_data_file  <- as.character(args[8])

library(bigsnpr)
library(tidyverse)

set.seed(2022)


################## Define paths ################ 

covariate_file = paste0("swapped_genotypes/covariates/covariates_", pheno.name, "_",
                        my_population ,"_sample_all_groups")

covariates_interactions = read.csv(covariate_file)
head(covariates_interactions)
names_covariates_interactions = colnames(covariates_interactions %>% dplyr::select(-eid))

# get last part after "_"
get_last_string <- function(string) {
  parts <- strsplit(string, "_")[[1]]
  return(tail(parts, 1))
}
filtered_col_names <- colnames(covariates_interactions)[!grepl("eid", colnames(covariates_interactions)) & !grepl("^rev", colnames(covariates_interactions))]
covar_names_pretty <- sapply(filtered_col_names, get_last_string)
covar_identifiyer = paste0(sort(covar_names_pretty), collapse = "_")

print("Covar identifier")
print(covar_identifiyer)

stats.path = paste0("UKB_stats/", my_population, 
                    "_res", resolution, "_", pheno.name, "_cov_", covar_identifiyer)

W.stats <- function(Z, knockoff) {
  importance <- abs(Z)
  z  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
  zk <- sum(importance[which(knockoff==TRUE)], na.rm=T)
  w <- z-zk
}


################ get real data ############

bed_bim_fam_filename <- real_data_file

# Attach the "bigSNP" object in R session
cat("Attaching bigSNP object... ")
obj.bigSNP.real.full <- snp_attach(paste0(bed_bim_fam_filename,".rds"))
cat("done.\n")

G.real = obj.bigSNP.real.full$genotypes

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

min_max_BP_by_group = map.real %>% 
  group_by(CHR, Group) %>% 
  summarize(min_BP = min(BP), max_BP = max(BP))

min_max_BP_by_group.file <- sprintf("%s_group_bp_min_max.txt", stats.path)
min_max_BP_by_group %>% write_delim(min_max_BP_by_group.file, delim=" ")
cat(sprintf("group info written on: %s\n", min_max_BP_by_group.file))


######### get phenotype ###########

cat("Reading phenotype file... ")
Phenotypes <- read.csv(pheno.file)
cat("done.\n")

ind.train <- which(!is.na(Phenotypes[[pheno.name]]))
y <- Phenotypes[[pheno.name]][ind.train]


# include all covariates, including those that we are looking for interactions for
baseline_covariate_names = c("age", "age_squared", "sex","PC1", "PC2", "PC3", "PC4", "PC5", "eid")
names_covariates_interactions = colnames(covariates_interactions %>% dplyr::select(-eid))

unique_baseline_covariate_names = baseline_covariate_names[!(baseline_covariate_names %in% names_covariates_interactions)] 
baseline_covariates <- Phenotypes %>% dplyr::select(all_of(unique_baseline_covariate_names))

covar.train = baseline_covariates %>% left_join(covariates_interactions, by = "eid") %>% dplyr::select(-eid)
covar.train = covar.train %>% dplyr::select(-starts_with("rev"))
covar.train = as.matrix(covar.train)[ind.train,,drop=F] 

######### run directly on real data ########### 

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
                          ind.train = ind.train,
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


selections_real <- Stats_real %>% knockoff.filter(fdr=fdr, offset=1) %>%
  arrange(CHR, Group) %>%
  arrange(desc(W))

print("Rows Stats real:")
print(nrow(Stats_real))

print(paste0("Number rejections without interactions directly with real data:", nrow(selections_real)))

########## running after pre-screening (on nonzero variants only) #######

# load nonzero snps from swapped lasso
coef.file <- sprintf("%s_coefficients_lasso_swapped.txt", stats.path)
coefs.swapped.lasso = read_delim(coef.file, delim=" ") %>% mutate(CHR_Group = paste0(CHR, "_", Group))
nonzero.chr.group = coefs.swapped.lasso %>%
  filter(Z != 0) %>%
  distinct(CHR_Group) %>% # any nonzero, does not matter whether .a or .b
  mutate(nonzero = 1) %>%
  select(CHR_Group, nonzero)

# save nonzero snps as bed/bim/fam
map.real = map.real %>%
  left_join(nonzero.chr.group, by = c("CHR_Group"))

index.nonzero.chr.group = which(map.real$nonzero == 1)
map.real.nonzero = map.real[index.nonzero.chr.group, ]

print("Dimensions nonzero map")
print(dim(map.real.nonzero))

# save nonzero SNPs
backingfile.name = paste0("temp/", my_population, "_res", resolution, "_" ,pheno.name, "_cov_", covar_identifiyer, "_nonzerorealsubset")

if(file.exists(paste0(backingfile.name, ".bk"))){
  file.remove(paste0(backingfile.name, ".bk"))
  file.remove(paste0(backingfile.name, ".rds"))
  file.remove(paste0(backingfile.name, ".bed"))
  file.remove(paste0(backingfile.name, ".bim"))
  file.remove(paste0(backingfile.name, ".fam"))
}

obj.bigSNP.real.subset.nonzero.file = snp_subset(obj.bigSNP.real.full, ind.col = index.nonzero.chr.group, backingfile = backingfile.name)
obj.bigSNP.real.subset.nonzero <- snp_attach(obj.bigSNP.real.subset.nonzero.file)
G.real.subset.nonzero = obj.bigSNP.real.subset.nonzero$genotypes

# Run lasso 
cat("Run lasso prescreened on real data ... \n")
scaler <- big_scale()
G.scale.screen <- scaler(G.real.subset.nonzero)
scaling.factors.screen <- G.scale.screen$scale

cat("Running lasso after screening... ")
lasso.fit.screen <- big_spLinReg(G.real.subset.nonzero,
                          y.train=y,
                          dfmax=dfmax, 
                          ncores=ncores,
                          ind.train = ind.train,
                          covar.train = covar.train, 
                          pf.covar=rep(0, ncol(covar.train)))

# Extract beta from each fold and combine them
beta.screen <- sapply(1:10, function(k) lasso.fit.screen[[1]][k][[1]]$beta)

# Separate the coefficients of the genetic variants from the coefficients of the covariates
total_number_coefficients.screen = nrow(beta.screen)
total_number_genotypes_in_model.screen = total_number_coefficients.screen - ncol(covar.train)
beta.variants.screen <- beta.screen[1:total_number_genotypes_in_model.screen,]

# Undo scaling of lasso coefficients
remaining_columns.screen = attr(lasso.fit.screen, "ind.col") # which genotypes were used for model fitting (index does not include covariates)
beta.variants.screen <- beta.variants.screen * scaling.factors.screen[remaining_columns.screen]
Beta_screened <- cbind(tibble(CHR=map.real.nonzero[remaining_columns.screen, ]$CHR,
                     SNP=map.real.nonzero[remaining_columns.screen, ]$SNP,
                     BP=map.real.nonzero[remaining_columns.screen, ]$BP,
                     Knockoff=map.real.nonzero[remaining_columns.screen, ]$Knockoff),
              as_tibble(beta.variants.screen)) %>% as_tibble()
colnames(Beta_screened) <- c("CHR", "SNP", "BP", "Knockoff",paste("K", seq(ncol(beta.variants.screen)),sep=""))
Beta_screened <- Beta_screened %>%
  mutate(Z=(K1+K2+K3+K4+K5+K6+K7+K8+K9+K10)/10,
         Nonzero=(K1!=0)+(K2!=0)+(K3!=0)+(K4!=0)+(K5!=0)+(K6!=0)+(K7!=0)+(K8!=0)+(K9!=0)+(K10!=0)) %>%
  select(CHR, SNP, BP, Knockoff,Z)

Stats_real_screened <- Beta_screened %>%
  inner_join(map.real.nonzero, by = c("CHR", "SNP", "BP", "Knockoff")) %>%
  select("CHR", "Group", "SNP", "BP", "Knockoff", "Z") %>%
  group_by(CHR, Group) %>%
  summarize(W = W.stats(abs(Z),Knockoff),
            Lead=which.max(abs(Z)), SNP.lead=SNP[Lead], BP.lead=BP[Lead],
            Size=n()) %>%
  ungroup() %>%
  select(CHR, Group, SNP.lead, BP.lead, Size, W)

# Run with real data only on subset
selections_real_screened <- Stats_real_screened %>% knockoff.filter(fdr=fdr, offset=1) %>%
  arrange(CHR, Group) %>%
  arrange(desc(W))

# Save list of discoveries
selections.out.file <- sprintf("%s_lasso_discoveries_no_interaction_real.txt", stats.path)
selections_real %>% write_delim(selections.out.file, delim=" ")
cat(sprintf("lasso discoveries real directly no interactions written on: %s\n", selections.out.file))

stats.out.file <- sprintf("%s_lasso_stats_no_interaction_real.txt", stats.path)
Stats_real %>% write_delim(stats.out.file, delim=" ")
cat(sprintf("lasso stats real directly no interactions written on: %s\n", stats.out.file))

selections.out.file <- sprintf("%s_lasso_discoveries_no_interaction_real_after_prescreen.txt", stats.path)
selections_real %>% write_delim(selections.out.file, delim=" ")
cat(sprintf("lasso discoveries real prescreened no interactions written on: %s\n", selections.out.file))

stats.out.file <- sprintf("%s_lasso_stats_no_interaction_real_after_prescreen.txt", stats.path)
Stats_real_screened %>% write_delim(stats.out.file, delim=" ")
cat(sprintf("lasso stats real prescreened no interactions written on: %s\n", stats.out.file))

