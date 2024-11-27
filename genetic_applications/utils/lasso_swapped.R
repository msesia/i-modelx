#!/usr/bin/env Rscript

# UKB GWAS
#
# Class: R script
#
# This script conducts the pre-screening 

args <- commandArgs(trailingOnly = TRUE)

my_population <- as.character(args[1])
resolution <- as.numeric(args[2])
pheno.name  <-  as.character(args[3])
dfmax <- as.numeric(args[4])
ncores <- as.numeric(args[5])
pheno.file <- as.character(args[6])

library(bigsnpr)
library(tidyverse)

print(version)

set.seed(2024)

######### load (ordered) swapped data ##############

ordered.bigsnp.filename = paste0("swapped_genotypes/chr_merged/UKB_swp_", my_population, 
                                 "_res", resolution, "_ordered_whr")

if(!file.exists(paste0(ordered.bigsnp.filename, ".bk"))) {
  ordered.bigsnp.filename.rds = snp_readBed2(paste0(ordered.bigsnp.filename, ".bed"))
}
cat("done.\n")

obj.bigSNP.swapped <- snp_attach(paste0(ordered.bigsnp.filename, ".rds"))
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

# Extract list of subjects
Subjects <- obj.bigSNP.swapped$fam %>% as_tibble()
colnames(Subjects) <- c("FID", "IID", "X1", "X2", "sex", "X3")
Subjects <- Subjects %>% select(FID, IID) %>%
  mutate(FID=as.character(FID), IID=as.character(IID))

print(dim(Subjects))

######## run lasso on swapped data #########

# Compute scaling factor for the genotypes
cat("Computing scaling factors for all variants... ")
scaler <- big_scale()
G.scale <- scaler(G.swapped)
scaling.factors <- G.scale$scale
cat("done.\n")

cat("Reading phenotype file... ")
Phenotypes <- read.csv(pheno.file)

print(head(Phenotypes))
print(paste0("Phenotype is: ", pheno.name))

# Extract response variable
sum(is.na(Phenotypes$whr))
ind.train <- which(!is.na(Phenotypes[[pheno.name]]))
cat(sprintf("%d individuals out of %d have missing phenotype.\n",
            nrow(Subjects)-length(ind.train), nrow(Subjects)))
y <- Phenotypes[[pheno.name]][ind.train]

# include all covariates, including those that we are looking for interactions for
baseline_covariate_names = c("age", "age_squared", "sex","PC1", "PC2", "PC3", "PC4", "PC5", "eid")

covariate_file=paste0("swapped_genotypes/covariates/covariates_", pheno.name, 
                      "_",my_population,"_sample_all_groups")

covariates_interactions = read.csv(covariate_file)
head(covariates_interactions)
names_covariates_interactions = colnames(covariates_interactions %>% dplyr::select(-eid))


# get covariates
get_last_string <- function(string) {
  parts <- strsplit(string, "_")[[1]]
  return(tail(parts, 1))
}
filtered_col_names <- colnames(covariates_interactions)[!grepl("eid", colnames(covariates_interactions)) & !grepl("^rev", colnames(covariates_interactions))]
covar_names_pretty <- sapply(filtered_col_names, get_last_string)
covar_identifiyer = paste0(sort(covar_names_pretty), collapse = "_")

print("Covar identifier")
print(covar_identifiyer)

unique_baseline_covariate_names = baseline_covariate_names[!(baseline_covariate_names %in% names_covariates_interactions)] 
baseline_covariates <- Phenotypes %>% dplyr::select(all_of(unique_baseline_covariate_names))

covar.train = baseline_covariates %>% left_join(covariates_interactions, by = "eid") %>% dplyr::select(-eid)
covar.train = covar.train %>% dplyr::select(-starts_with("rev"))
covar.train = as.matrix(covar.train)[ind.train,,drop=F] 

# Make sure that the rows of the genotypes match the rows of the phenotypes
# all are unrelated
if(sum(Phenotypes$eid != Subjects$FID) > 0) {
  stop("Subjects not in correct order!")
}

cat(sprintf("Fitting sparse linear regression with %d observations, %d variants and %d covariates... ",
            length(y), ncol(G.swapped), ncol(covar.train)))
lasso.fit <- big_spLinReg(G.swapped, y.train=y, ind.train=ind.train, covar.train=covar.train,
                          dfmax=dfmax, ncores=ncores) #, pf.covar=rep(0, ncol(covar.train))


print(paste0("dfmax is ", dfmax))

# Extract beta from each fold and combine them
cat("Extracting regression coefficients... ")
beta <- sapply(1:10, function(k) lasso.fit[[1]][k][[1]]$beta)

# Separate the coefficients of the genetic variants from the coefficients of the covariates
beta.variants <- beta[1:ncol(G.swapped),]
beta.covariates <- beta[(ncol(G.swapped)+1):nrow(beta),]

# Undo scaling of lasso coefficients
beta.variants <- beta.variants * scaling.factors
Beta <- cbind(tibble(CHR=map.swapped$CHR,
                     SNP=map.swapped$SNP, BP=map.swapped$physical.pos, knock_ind = map.swapped$knock_ind),
              as_tibble(beta.variants)) %>% as_tibble()
colnames(Beta) <- c("CHR", "SNP", "BP", "knock_ind",paste("K", seq(ncol(beta.variants)),sep=""))
Beta <- Beta %>%
  mutate(Z=(K1+K2+K3+K4+K5+K6+K7+K8+K9+K10)/10,
         Nonzero=(K1!=0)+(K2!=0)+(K3!=0)+(K4!=0)+(K5!=0)+(K6!=0)+(K7!=0)+(K8!=0)+(K9!=0)+(K10!=0))

print(paste0("Number nonzero beta: ", sum(Beta$Z != 0)))

Beta = Beta %>%
  select(CHR, SNP, BP, knock_ind, Z)

Beta <- Beta %>% 
  left_join(map.swapped, by = c("CHR", "SNP", "knock_ind"))

######## select any nonzero snps ###########

#(regardless of whether (swapped) knockoff or original)
nonzero.chr.group = Beta %>% 
  filter(Z != 0) %>% 
  distinct(CHR_Group) %>% # any nonzero, does not matter whether .a or .b
  mutate(nonzero = 1) %>% 
  select(CHR_Group, nonzero)

print(paste0("Number nonzero CHR Group: ", nrow(nonzero.chr.group)))

nonzero.snp = Beta %>% 
  filter(Z != 0) %>% 
  distinct(SNP) %>% # any nonzero, does not matter whether .a or .b
  mutate(nonzero = 1) %>% 
  select(SNP, nonzero)

print(paste0("Number nonzero SNPS: ", nrow(nonzero.snp)))

# save nonzero snps / group as bed/bim/fam
map.swapped = map.swapped %>% 
  left_join(nonzero.snp, by = c("SNP"))

index.nonzero.snp = which(map.swapped$nonzero == 1)

# save nonzero SNPs
backingfile.name = paste0("swapped_genotypes/chr_merged/prescreened/UKB_swp_", 
                          my_population, "_res", resolution, "_" ,pheno.name, "_cov_", covar_identifiyer, "_ordered", "_prescreened")

if(file.exists(paste0(backingfile.name, ".bk"))){
  file.remove(paste0(backingfile.name, ".bk"))
  file.remove(paste0(backingfile.name, ".rds"))
  file.remove(paste0(backingfile.name, ".bed"))
  file.remove(paste0(backingfile.name, ".bim"))
  file.remove(paste0(backingfile.name, ".fam"))
}

obj.bigSNP.swapped.prescreened.file = snp_subset(obj.bigSNP.swapped, ind.col = index.nonzero.snp, 
                                                 backingfile = backingfile.name)

obj.bigSNP.swapped.prescreened <- snp_attach(obj.bigSNP.swapped.prescreened.file)


# save as bed/bim/fam 
output <- structure(list(genotypes = obj.bigSNP.swapped.prescreened$genotypes, 
                         fam = obj.bigSNP.swapped.prescreened$fam,
                         map = obj.bigSNP.swapped.prescreened$map),
                    class = "bigSNP")

snp_writeBed(output, paste0(backingfile.name, ".bed"))

print(paste0("Number nonzero screened genotypes: ", ncol(obj.bigSNP.swapped.prescreened$genotypes)))

# save betas
stats.path = paste0("UKB_stats/", my_population, 
                    "_res", resolution, "_", pheno.name, "_cov_", covar_identifiyer)
coef.file <- sprintf("%s_coefficients_lasso_swapped.txt", stats.path)
Beta %>% write_delim(coef.file, delim=" ")
cat(sprintf("Betas written on:\n%s\n", coef.file))