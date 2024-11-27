#!/usr/bin/env Rscript

# UK Biobank GWAS
#
# Class: R script
#
# Run lasso and obtain knockoff statistics within a specified partition

# Input arguments
args <- commandArgs(trailingOnly=TRUE)
my_population     <- as.character(args[1])
resolution       <- as.integer(args[2])
pheno.name       <- as.character(args[3])
interaction_character <- as.character(args[4])
environment       <- as.integer(args[5])
pheno.file       <- as.character(args[6])
dfmax       <- as.integer(args[7])
ncores       <- as.integer(args[8])
suffix       <- as.character(args[9])
num.int       <- as.integer(args[10])
real_data_file <- as.character(args[11])

set.seed(2022)

library(bigsnpr)
library(tidyverse)

print(version)

print(paste0("Population: ", my_population))
print(paste0("Resolution: ", resolution))
print(paste0("Phenotype: ", pheno.name))
print(paste0("Interaction: ", interaction_character))
print(paste0("Environment: ", environment))
print(paste0("Pheno file: ", pheno.file))
print(paste0("dfmax: ", dfmax))
print(paste0("ncores: ", ncores))

W.stats <- function(Z, knockoff) {
  importance <- abs(Z)
  z  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
  zk <- sum(importance[which(knockoff==TRUE)], na.rm=T)
  w <- z-zk
}

# set path where to save results
stats.path = paste0("UKB_stats/lasso_env_stats")

# get covariates
covariate_file = paste0("covariates/covariates_", pheno.name, "_",
                        my_population ,"_sample_all_groups")
covariates_interactions = read.csv(covariate_file)

# get last part after "_" for covariate identification
get_last_string <- function(string) {
  parts <- strsplit(string, "_")[[1]]
  return(tail(parts, 1))
}

filtered_col_names <- colnames(covariates_interactions)[!grepl("eid", colnames(covariates_interactions)) & !grepl("^rev", colnames(covariates_interactions))]
covar_names_pretty <- sapply(filtered_col_names, get_last_string)
covar_identifiyer = paste0(sort(covar_names_pretty), collapse = "_")

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

######## get swapped data #############

swapped_genotypes_path <- "swapped_genotypes"

swapped.genotypes.filename <-  paste0(swapped_genotypes_path, "/chr_merged/UKB_swp_", my_population, 
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


# Extract list of subjects
Subjects <- obj.bigSNP.swapped$fam %>% as_tibble()
colnames(Subjects) <- c("FID", "IID", "X1", "X2", "sex", "X3")
Subjects <- Subjects %>% select(FID, IID) %>%
  mutate(FID=as.character(FID), IID=as.character(IID))


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

########### get environment and chr group info ##################

environment_chr_group_info = readRDS(paste0("lasso_each_env/",
                                            my_population, "_res", resolution, "_", pheno.name, "_cov_", covar_identifiyer, "_", 
                                            interaction_character, "_env", environment, "_", suffix, "_", num.int))

relevant_chr_group_index = environment_chr_group_info$relevant_chr_group_index
environment_index_obs = environment_chr_group_info$environment_index_obs

bk.filename = paste0(swapped_genotypes_path, "tmp/lasso_each_env_",
                     my_population, "_res", resolution, "_", pheno.name, "_cov_", covar_identifiyer, 
                     "_", interaction_character,"_env" ,environment, "_", suffix, "_", num.int)
if(file.exists(paste0(bk.filename, "_swapped.bk"))) {
  file.remove(paste0(bk.filename, "_swapped.bk"))
}
if(file.exists(paste0(bk.filename, "_real.bk"))) {
  file.remove(paste0(bk.filename, "_real.bk"))
}

# get chr-group from chr-group index
relevant_chr_group = unique(map.real[relevant_chr_group_index, ]$CHR_Group)

print("Number of unique chr groups")
print(length(unique(relevant_chr_group)))

######### get phenotype ###########

cat("Reading phenotype file... ")
Phenotypes <- read.csv(pheno.file)
cat("done.\n")

# Extract response variable
ind.train <- which(!is.na(Phenotypes[[pheno.name]]))
cat(sprintf("%d individuals out of %d have missing phenotype.\n",
            nrow(Subjects)-length(ind.train), nrow(Subjects)))

y <- Phenotypes[[pheno.name]][intersect(environment_index_obs, ind.train)]

if(sum(is.na(y)) > 0) {
  stop("Y has missing values")
}

baseline_covariate_names = c("age", "age_squared", "sex","PC1", "PC2", "PC3", "PC4", "PC5", "eid")
names_covariates_interactions = colnames(covariates_interactions %>% dplyr::select(-eid))

unique_baseline_covariate_names = baseline_covariate_names[!(baseline_covariate_names %in% names_covariates_interactions)] 
baseline_covariates <- Phenotypes %>% dplyr::select(all_of(unique_baseline_covariate_names))

covar.train = baseline_covariates %>% left_join(covariates_interactions, by = "eid") %>% dplyr::select(-eid)
covar.train = covar.train %>% dplyr::select(-starts_with("rev"))
covar.train = covar.train[intersect(environment_index_obs, ind.train),,drop=F] 

# filter out covariates that do not have any variance
variance_each_covar = sapply(covar.train, function(x) c(var=var(x)))
print("null variance covars ")
variance_each_covar[which(variance_each_covar == 0)]
names(variance_each_covar) = NULL

if(sum(variance_each_covar == 0) > 0) {
  covar.train = as.matrix(covar.train[, -which(variance_each_covar == 0)])
} else {
  covar.train = as.matrix(covar.train)
}

head(covar.train)
print(dim(covar.train))


########### subset genotype, phenotype files ##########

# subset genotype files to only those with observations in the current environment

cat("Subsetting real and swapped data ... ")

G.swapped.subset.env = big_copy(G.swapped, ind.row = intersect(environment_index_obs, ind.train), backingfile = paste0(bk.filename, "_swapped"))
G.real.subset.env = big_copy(G.real, ind.row = intersect(environment_index_obs, ind.train), backingfile = paste0(bk.filename, "_real"))
cat("done.\n")

# replace each relevant chr-group with its real counterpart
cat("Replacing swapped with real data for relevant chr groups ... ")
big_apply(G.swapped.subset.env, function(G.swapped.subset.env, G.real.subset.env, ind) {
  
  G.swapped.subset.env[, ind] <- G.real.subset.env[, ind]
  
  # return nothing
  NULL
  
}, block.size = 1, a.combine = 'cbind',
ind = relevant_chr_group_index,
G.real.subset.env = G.real.subset.env)
cat("done.\n")

###### run lasso #########
scaler <- big_scale()
G.scale <- scaler(G.swapped.subset.env)
scaling.factors <- G.scale$scale

# Run lasso
cat(sprintf("Fitting sparse linear regression with %d observations, %d variants and %d covariates... ",
            length(y), ncol(G.swapped.subset.env), ncol(covar.train)))
lasso.fit <- big_spLinReg(G.swapped.subset.env, y.train=y, dfmax=dfmax,
                          ncores=ncores, covar.train =covar.train) # , covar.train = covar.train , pf.covar=rep(0, ncol(covar.train))

# Extract beta from each fold and combine them
beta <- sapply(1:10, function(k) lasso.fit[[1]][k][[1]]$beta)

# Separate the coefficients of the genetic variants from the coefficients of the covariates
total_number_coefficients = nrow(beta)
total_number_genotypes_in_model = total_number_coefficients - ncol(covar.train)
beta.variants <- beta[1:total_number_genotypes_in_model,]

# Undo scaling of lasso coefficients
remaining_columns = attr(lasso.fit, "ind.col") # which genotypes were used for model fitting (index does not include covariates)
beta.variants <- beta.variants * scaling.factors[remaining_columns]
Beta <- cbind(tibble(CHR=map.swapped[remaining_columns, ]$CHR,
                     SNP=map.swapped[remaining_columns, ]$SNP,
                     physical.pos=map.swapped[remaining_columns, ]$physical.pos,
                     Knockoff=map.swapped[remaining_columns, ]$Knockoff),
              as_tibble(beta.variants)) %>% as_tibble()
colnames(Beta) <- c("CHR", "SNP", "physical.pos", "Knockoff",paste("K", seq(ncol(beta.variants)),sep=""))
Beta <- Beta %>%
  mutate(Z=(K1+K2+K3+K4+K5+K6+K7+K8+K9+K10)/10,
         Nonzero=(K1!=0)+(K2!=0)+(K3!=0)+(K4!=0)+(K5!=0)+(K6!=0)+(K7!=0)+(K8!=0)+(K9!=0)+(K10!=0)) %>%
  select(CHR, SNP, physical.pos, Knockoff,Z)

# Extract the estimated coefficients
Stats <- Beta %>%
  inner_join(map.swapped, by = c("CHR", "SNP", "physical.pos", "Knockoff")) %>%
  select("CHR", "Group", "SNP", "physical.pos", "Knockoff", "Z") %>%
  group_by(CHR, Group) %>%
  summarize(W = W.stats(abs(Z),Knockoff),
            Lead=which.max(abs(Z)), SNP.lead=SNP[Lead], BP.lead=physical.pos[Lead],
            Size=n()) %>%
  ungroup() %>%
  select(CHR, Group, SNP.lead, BP.lead, Size, W) %>%
  mutate(Interaction = interaction_character,
         environment = environment)

Stats_filtered = Stats %>%
  mutate(CHR_Group = paste0(CHR, "_", Group)) %>%
  filter(CHR_Group %in% relevant_chr_group)
cat("done.\n")


Beta_filtered = Beta %>%
  inner_join(map.swapped, by = c("CHR", "SNP", "physical.pos", "Knockoff")) %>%
  select("CHR", "Group", "SNP", "physical.pos", "Knockoff", "Z") %>%
  mutate(CHR_Group = paste0(CHR, "_", Group)) %>%
  filter(CHR_Group %in% relevant_chr_group)


print("Dimensions stats filtered")
print(dim(Stats_filtered))


print("Dimensions stats")
print(dim(Stats))

# save stats
Stats_filtered_file <-   paste0(stats.path, "_",my_population, "_res", resolution, "_", pheno.name, "_cov_", covar_identifiyer, "_", 
                                interaction_character, "_env", environment, "_", suffix, "_", num.int ,".txt")
Stats_filtered %>% write_delim(Stats_filtered_file, delim=" ")
cat(sprintf("Statistics written on: %s\n", Stats_filtered_file))

beta_filtered_file <-   paste0(stats.path, "_",my_population, "_res", resolution, "_", pheno.name, "_cov_", covar_identifiyer, "_", 
                               interaction_character, "_env", environment, "_", suffix, "_", num.int , "_beta_filtered" ,".txt")
Beta_filtered %>% write_delim(beta_filtered_file, delim=" ")
cat(sprintf("Betas filtered written on: %s\n", beta_filtered_file))

Stats_file <-   paste0(stats.path, "_",my_population, "_res", resolution, "_", pheno.name, "_cov_", covar_identifiyer, "_", 
                       interaction_character, "_env", environment, "_", suffix, "_", num.int ,"_unfiltered.txt")
Stats %>% write_delim(Stats_file, delim=" ")
cat(sprintf("Statistics unfiltered written on: %s\n", Stats_file))

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

# knockoff filter only on partition
fdr = 0.1
selections <- Stats_filtered %>% knockoff.filter(fdr=fdr, offset=1) %>% 
  arrange(desc(W)) 

print("Number rejections this environment only")
print(dim(selections))

selections_file <-   paste0(stats.path, "_",my_population, "_res", resolution, "_", pheno.name, "_cov_", covar_identifiyer, "_", 
                            interaction_character, "_env", environment, "_", suffix, "_", num.int ,"_selections.txt")
selections %>% write_delim(selections_file, delim=" ")
cat(sprintf("selections filtered written on: %s\n", selections_file))



