#!/usr/bin/env Rscript

# Simulation based on UK Biobank genotypes
#
# Class: R script
#
# Simulate performance of fixed-LKF (1 covariate, i.e. two subgroups)


# Input arguments
args <- commandArgs(trailingOnly=TRUE)
sparsity <- as.numeric(args[1])
amp <- as.numeric(args[2])
fdr  <- as.numeric(args[3])
ncores <- as.numeric(args[4])
dfmax  <- as.numeric(args[5])
dfmax_lasso_all  <- as.numeric(args[6])
B  <- as.numeric(args[7])
propmain <- as.numeric(args[8])
prop_z0 <- as.numeric(args[9])
prop_z1 <- as.numeric(args[10])
random_sample_prop <- as.numeric(args[11])


library(bigsnpr)
library(tidyverse)

print(version)

my_population = "whitenonbritish"
resolution = 1

num_covariates = 1 # number of covariates

fdp = c()
fdp_interactions = c()
power_cond_association = c()
homogeneity = c()
homogeneity_complement = c()

fdp_interactions = c()
power_interactions = c()
fdp_main = c()
power_main = c()

fdp_z0 = c()
fdp_z1 = c()

power_z0 = c()
power_z1 = c()

power_cond_association_deconstruct_int = c()
power_cond_association_deconstruct_main = c()
power_cond_association_deconstruct_main_z1 = c()
power_cond_association_deconstruct_int_z1 = c()
power_cond_association_deconstruct_main_z0 = c()
power_cond_association_deconstruct_int_z0 = c()

homogeneity_z0 = c()
homogeneity_z1 = c()

heritability = c()

print(paste0("Number of iterations: ", B))

# load real genotypes
bed_bim_fam_filename <- paste0("ukb_", 
                               my_population, 
                               "_plh_gen_merged_filtered_res", 
                               resolution)


# Attach the "bigSNP" object in R session
cat("Attaching bigSNP object... ")
obj.bigSNP.real.full <- snp_attach(paste0(bed_bim_fam_filename,".rds"))
cat("done.\n")

G.real = obj.bigSNP.real.full$genotypes

# get list of variants 
map.real.main = obj.bigSNP.real.full$map %>% as_tibble()
colnames(map.real.main) <- c("CHR", "SNP", "gd", "BP", "a1", "a2")
map.real.main <- map.real.main %>% dplyr::select(CHR, SNP, BP) %>% 
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

map.real.main = map.real.main %>% left_join(grp.file, by = c("CHR", "SNP")) %>% mutate(CHR_Group = paste0(CHR, "_", Group))


# get real genotypes and map for chr 21 only
map.real.chr21.main = map.real.main %>% filter(CHR == 21)

G.real.chr21 = big_copy(G.real, ind.col = which(map.real.main$CHR == 21),
                        backingfile = tempfile(tmpdir = "tmp"))

# subset to chr 21 only 
subset.index.chr21.orig =  which(map.real.main$CHR == 21 & map.real.main$Knockoff == FALSE)
map.real.chr21.orig.main = map.real.main[subset.index.chr21.orig, ]
G.real.chr21.orig = big_copy(G.real, ind.col = subset.index.chr21.orig,
                             backingfile = tempfile(tmpdir = "tmp"))

# get random sample
if(!is.na(random_sample_prop)) {
  set.seed(123)
  ind.sample = sample(1:nrow(G.real.chr21), round(random_sample_prop*nrow(G.real.chr21))) # prop is with respect to chr 21 SNPs only
  
  G.real.chr21 = big_copy(G.real.chr21, ind.row = ind.sample,
                          backingfile = tempfile(tmpdir = "tmp"))
  
  
  G.real.chr21.orig = big_copy(G.real.chr21.orig, ind.row = ind.sample,
                               backingfile = tempfile(tmpdir = "tmp"))
  
}



for(b in 1:B) {
  
  print(paste0("iter ", b))
  
  set.seed(2*b + 1)
  
  ############ simulate y based on real genotypes (chr 21, wnb) ######################
  
  p = ncol(G.real.chr21.orig) 
  n = nrow(G.real.chr21.orig)
  unique_chr_group = sort(unique(map.real.chr21.orig.main$CHR_Group))
  number_unique_chr_group = length(unique(map.real.chr21.orig.main$CHR_Group))
  
  
  map.real.chr21.orig = map.real.chr21.orig.main
  map.real = map.real.main
  map.real.chr21 = map.real.chr21.main
  
  ############## simulate interactions ################## 
  
  # sample covariates
  covariates = c()
  for(q in 1:num_covariates) {
    Z = rbinom(n, 1, 1/2)
    covariates = cbind(covariates, Z)
  }
  colnames(covariates) = paste0("Z", seq(1, num_covariates))
  covariates_z11 = covariates[, 1]
  covariates_z10 = 1 - covariates[, 1]
  
  if(num_covariates == 2) {
    covariates_z21 = covariates[, 2]
    covariates_z20 = 1 - covariates[, 2]
  }
  
  
  # sample nonzero chr-group (sparsity on snp level is approx sparsity / 4) 
  number_nonzero_chr_group = round(sparsity*number_unique_chr_group)
  nonzero_chr_groups = sample(unique_chr_group, number_nonzero_chr_group)
  
  # from each nonzero chr group, randomly sample a single snp
  sampled_nonzero_snps <- map.real.chr21.orig %>%
    filter(CHR_Group %in% nonzero_chr_groups) %>%
    group_by(CHR_Group) %>%
    slice_sample(n = 1) %>% 
    select(CHR_Group, SNP) %>% 
    mutate(nonzero = 1)
  
  # randomly assign each non-zero snps to interaction or main effects
  p_z0 = prop_z0*(1 - propmain)
  p_z1 = prop_z1*(1 - propmain)
  
  sample_nonzero_interaction_vector = sample(1:3, number_nonzero_chr_group, replace = TRUE, prob = c(propmain, p_z1, p_z0))

  if(propmain == 0 & (p_z1 != 1 & p_z0 != 1)) {
    sample_nonzero_interaction = model.matrix(~0+as.factor(sample_nonzero_interaction_vector)) #one-hot-encode
    sample_nonzero_interaction = cbind(rep(0, number_nonzero_chr_group), sample_nonzero_interaction)
    colnames(sample_nonzero_interaction) = c("nonzero_all","nonzero_z11", "nonzero_z10")
  } 
  
  if(propmain == 0 & (p_z1 == 1 | p_z0 == 1)) {
    if(p_z1 == 1) {
      sample_nonzero_interaction = cbind(rep(0, number_nonzero_chr_group), 
                                         rep(1, number_nonzero_chr_group), 
                                         rep(0, number_nonzero_chr_group))
    } 
    
    
    if(p_z0 == 1) {
      sample_nonzero_interaction = cbind(rep(0, number_nonzero_chr_group), 
                                         rep(0, number_nonzero_chr_group), 
                                         rep(1, number_nonzero_chr_group))
    }
    colnames(sample_nonzero_interaction) = c("nonzero_all","nonzero_z11", "nonzero_z10")
  } 
  
  if(propmain == 1) {
    sample_nonzero_interaction = cbind(sample_nonzero_interaction_vector, 
                                       rep(0, number_nonzero_chr_group), 
                                       rep(0, number_nonzero_chr_group))
    colnames(sample_nonzero_interaction) = c("nonzero_all","nonzero_z11", "nonzero_z10")
  } 
  
  if(propmain != 0 &  propmain != 1) {
    if(length(unique(sample_nonzero_interaction_vector)) < 3) {
      missing_index = setdiff(seq(1:3), unique(sample_nonzero_interaction_vector))
      
      if(length(missing_index) == 1) {
        sample_nonzero_interaction_matrix_non_missing = model.matrix(~0+as.factor(sample_nonzero_interaction_vector)) #one-hot-encode
        
        sample_nonzero_interaction = matrix(0, nrow = nrow(sample_nonzero_interaction_matrix_non_missing), ncol= 3)
        sample_nonzero_interaction[, seq(1:3)[-missing_index]] <- sample_nonzero_interaction_matrix_non_missing
        colnames(sample_nonzero_interaction) = c("nonzero_all","nonzero_z11", "nonzero_z10")      
      } else {
        sample_nonzero_interaction = matrix(0, nrow = length(sample_nonzero_interaction_vector), ncol= 3)
        sample_nonzero_interaction[, seq(1:3)[-c(missing_index)]] <- sample_nonzero_interaction_vector
        colnames(sample_nonzero_interaction) = c("nonzero_all","nonzero_z11", "nonzero_z10")     
      }
      
      
      
    } else {
      sample_nonzero_interaction = model.matrix(~0+as.factor(sample_nonzero_interaction_vector)) #one-hot-encode
      colnames(sample_nonzero_interaction) = c("nonzero_all","nonzero_z11", "nonzero_z10")
    }
  }
  
  colSums(sample_nonzero_interaction)
  
  
  sampled_nonzero_snps = cbind(sampled_nonzero_snps, sample_nonzero_interaction)
  
  sampled_nonzero_snps = sampled_nonzero_snps %>% 
    mutate(Interaction = ifelse(nonzero_all == 1, "none", "Z1"), 
           environment = ifelse(nonzero_z11 == 1, 2, 1))
  
  print("number nonzero by interaction")
  colSums(sampled_nonzero_snps[, 3:6])
  
  map.real.chr21.orig = map.real.chr21.orig %>% left_join(sampled_nonzero_snps, by = c("CHR_Group", "SNP"))
  index_nonzero = which(map.real.chr21.orig$nonzero == 1)
  index_nonzero_all = which(map.real.chr21.orig$nonzero_all == 1)
  index_nonzero_int_z11 = which(map.real.chr21.orig$nonzero_z11 == 1)
  index_nonzero_int_z10 = which(map.real.chr21.orig$nonzero_z10 == 1)
  
  # create coefficients
  factor_all = nrow(covariates) / nrow(covariates)
  factor_z11 = nrow(covariates) / sum(covariates[, 1] == 1)
  factor_z10 = nrow(covariates) / sum(covariates[, 1] == 0)

  beta_all <- amp * sqrt(factor_all) * (1:p %in% index_nonzero_all)*sign(rnorm(p,0,1))/sqrt(n)
  beta_z11 <- amp * sqrt(factor_z11) * (1:p %in% index_nonzero_int_z11)*sign(rnorm(p,0,1))/sqrt(n)
  beta_z10 <- amp * sqrt(factor_z10) * (1:p %in% index_nonzero_int_z10)*sign(rnorm(p,0,1))/sqrt(n)
  
  # generate y
  y_all <- big_prodMat(G.real.chr21.orig, as.matrix(beta_all)) 
  y_intz11 <- big_prodMat(G.real.chr21.orig, as.matrix(beta_z11)) * covariates_z11
  y_intz10 <- big_prodMat(G.real.chr21.orig, as.matrix(beta_z10)) * covariates_z10
  
  x_beta = y_all + y_intz11 + y_intz10 
  y = x_beta + rnorm(n)
  
  heritability[b] = var(x_beta) / var(y)
  
  # effect matrix
  ite = matrix(0, nrow = n, ncol = p)
  for(j in 1:p) {
    ite[,j ] <- beta_all[j]*rep(1, n) + beta_z10[j]*covariates_z10 +  beta_z11[j]*covariates_z11
  }
  
  # save simulated y, clear entire environment
  rm(list=setdiff(ls(), c("y", "covariates",  "map.real.chr21.orig.main", "map.real.chr21.main", "map.real.chr21", "map.real.chr21.orig",  "map.real",  "map.real.main", 
                          "G.real.chr21", "G.real.chr21.orig", "my_population", "sparsity", "fdr", "amp", 
                          "ncores", "dfmax", "B", "b", "resolution",  "num_covariates", "power_cond_association_deconstruct_int", "power_cond_association_deconstruct_main",
                           "dfmax_lasso_all",  "sample_nonzero_interaction", "sampled_nonzero_snps",
                          "fdp", "power_cond_association", "propmain", "ite", "homogeneity", 
                          "power_main", "fdp_main", "power_interactions", "fdp_interactions", "heritability", 
                          "power_cond_association_deconstruct_main_z0", "power_cond_association_deconstruct_int_z0",
                          "power_cond_association_deconstruct_main_z1", "power_cond_association_deconstruct_int_z1",
                          "homogeneity_z0", "homogeneity_z1",
                          "fdp_z0", "fdp_z1", "power_z0", "power_z1", "prop_z1", "prop_z0", "random_sample_prop")))
  
  # calculate sum of importance statistics from swapped lasso
  W.stats <- function(Z, knockoff) {
    importance <- abs(Z)
    z  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
    zk <- sum(importance[which(knockoff==TRUE)], na.rm=T)
    w <- z-zk
  }
  
  
  ######################### Subset Z = 0 ##########################################
  
  
  ind.z0 = which(covariates[, 1] == 0)
  
  G.real.chr21.z0 = big_copy(G.real.chr21, ind.row = ind.z0)
  y.z0 = y[ind.z0]
  
  covariates.z0 = covariates[ind.z0, , drop = FALSE]
  
  # Compute scaling factor for the genotypes
  cat("Computing scaling factors for all variants... ")
  scaler <- big_scale()
  G.scale <- scaler(G.real.chr21.z0)
  scaling.factors.z0 <- G.scale$scale
  cat("done.\n")
  
  # Run lasso on swapped genotypes for pre-screening 
  cat(sprintf("Fitting sparse linear regression with %d observations, %d variants and NO covariates... ",
              length(y.z0), ncol(G.real.chr21.z0)))
  lasso.fit.z0 <- big_spLinReg(G.real.chr21.z0, y.train=y.z0, dfmax=dfmax, ncores=ncores)
  
  # Extract beta from each fold and combine them
  cat("Extracting regression coefficients... ")
  beta.z0 <- sapply(1:10, function(k) lasso.fit.z0[[1]][k][[1]]$beta)
  
  # get coefficients
  total_number_genotypes_in_model.z0 = nrow(beta.z0) 
  beta.variants.z0 <- beta.z0[1:total_number_genotypes_in_model.z0,]
  
  # Undo scaling of lasso coefficients
  remaining_columns.z0 = attr(lasso.fit.z0, "ind.col") # which genotypes were used for model fitting (index does not include covariates)
  beta.variants.z0 <- beta.variants.z0 * scaling.factors.z0[remaining_columns.z0]
  
  Beta.z0 <- cbind(tibble(CHR=map.real.chr21[remaining_columns.z0,]$CHR,
                       SNP=map.real.chr21[remaining_columns.z0,]$SNP, BP=map.real.chr21[remaining_columns.z0,]$BP, 
                       Knockorr = map.real.chr21[remaining_columns.z0,]$Knockoff),
                as_tibble(beta.variants.z0)) %>% as_tibble()
  colnames(Beta.z0) <- c("CHR", "SNP", "BP", "Knockoff",paste("K", seq(ncol(beta.variants.z0)),sep=""))
  Beta.z0 <- Beta.z0 %>%
    mutate(Z=(K1+K2+K3+K4+K5+K6+K7+K8+K9+K10)/10,
           Nonzero=(K1!=0)+(K2!=0)+(K3!=0)+(K4!=0)+(K5!=0)+(K6!=0)+(K7!=0)+(K8!=0)+(K9!=0)+(K10!=0)) %>%
    select(CHR, SNP, BP, Knockoff, Z)
  
  Beta.z0 <- Beta.z0 %>% 
    left_join(map.real.chr21, by = c("CHR", "SNP", "Knockoff", "BP"))
  
  Stats.z0 <- Beta.z0 %>%
    select("CHR", "Group", "CHR_Group", "SNP", "BP", "Knockoff", "Z") %>%
    group_by(CHR, Group) %>%
    summarize(W = W.stats(abs(Z),Knockoff),
              Lead=which.max(abs(Z)), SNP.lead=SNP[Lead], BP.lead=BP[Lead],
              Size=n()) %>%
    ungroup() %>%
    select(CHR, Group, SNP.lead, BP.lead, Size, W) %>%
    mutate(CHR_Group = paste0(CHR, "_", Group))  %>% 
    mutate(Interaction.selected = "Z1", 
           environment.selected = 1)
  
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
  
  
  # lasso results with interaction
  selections.lasso.z0 <- Stats.z0 %>% knockoff.filter(fdr=fdr, offset=1) %>% 
    arrange(CHR, Group) %>% 
    arrange(desc(W)) 
  
  print(paste0("Number selections lasso z0: ", nrow(selections.lasso.z0)))

  
  ######################### Subset Z = 1 ##########################################
  
  
  ind.z1 = which(covariates[, 1] == 1)
  
  G.real.chr21.z1 = big_copy(G.real.chr21, ind.row = ind.z1)
  y.z1 = y[ind.z1]
  
  covariates.z1 = covariates[ind.z1, , drop = FALSE]
  
  # Compute scaling factor for the genotypes
  cat("Computing scaling factors for all variants... ")
  scaler <- big_scale()
  G.scale <- scaler(G.real.chr21.z1)
  scaling.factors.z1 <- G.scale$scale
  cat("done.\n")
  
  # Run lasso on swapped genotypes for pre-screening 
  cat(sprintf("Fitting sparse linear regression with %d observations, %d variants and NO covariates... ",
              length(y.z1), ncol(G.real.chr21.z1)))
  lasso.fit.z1 <- big_spLinReg(G.real.chr21.z1, y.train=y.z1, dfmax=dfmax, ncores=ncores)
  
  # Extract beta from each fold and combine them
  cat("Extracting regression coefficients... ")
  beta.z1 <- sapply(1:10, function(k) lasso.fit.z1[[1]][k][[1]]$beta)
  
  # get coefficients
  total_number_genotypes_in_model.z1 = nrow(beta.z1) 
  beta.variants.z1 <- beta.z1[1:total_number_genotypes_in_model.z1,]
  
  # Undo scaling of lasso coefficients
  remaining_columns.z1 = attr(lasso.fit.z1, "ind.col") # which genotypes were used for model fitting (index does not include covariates)
  beta.variants.z1 <- beta.variants.z1 * scaling.factors.z1[remaining_columns.z1]
  
  Beta.z1 <- cbind(tibble(CHR=map.real.chr21[remaining_columns.z1,]$CHR,
                          SNP=map.real.chr21[remaining_columns.z1,]$SNP, BP=map.real.chr21[remaining_columns.z1,]$BP, 
                          Knockorr = map.real.chr21[remaining_columns.z1,]$Knockoff),
                   as_tibble(beta.variants.z1)) %>% as_tibble()
  colnames(Beta.z1) <- c("CHR", "SNP", "BP", "Knockoff",paste("K", seq(ncol(beta.variants.z1)),sep=""))
  Beta.z1 <- Beta.z1 %>%
    mutate(Z=(K1+K2+K3+K4+K5+K6+K7+K8+K9+K10)/10,
           Nonzero=(K1!=0)+(K2!=0)+(K3!=0)+(K4!=0)+(K5!=0)+(K6!=0)+(K7!=0)+(K8!=0)+(K9!=0)+(K10!=0)) %>%
    select(CHR, SNP, BP, Knockoff, Z)
  
  Beta.z1 <- Beta.z1 %>% 
    left_join(map.real.chr21, by = c("CHR", "SNP", "Knockoff", "BP"))
  
  Stats.z1 <- Beta.z1 %>%
    select("CHR", "Group", "CHR_Group", "SNP", "BP", "Knockoff", "Z") %>%
    group_by(CHR, Group) %>%
    summarize(W = W.stats(abs(Z),Knockoff),
              Lead=which.max(abs(Z)), SNP.lead=SNP[Lead], BP.lead=BP[Lead],
              Size=n()) %>%
    ungroup() %>%
    select(CHR, Group, SNP.lead, BP.lead, Size, W) %>%
    mutate(CHR_Group = paste0(CHR, "_", Group))  %>% 
    mutate(Interaction.selected = "Z1", 
           environment.selected = 2)
  
  
  # lasso results with interaction
  selections.lasso.z1 <- Stats.z1 %>% knockoff.filter(fdr=fdr, offset=1) %>% 
    arrange(CHR, Group) %>% 
    arrange(desc(W))  
  print(paste0("Number selections lasso z1: ", nrow(selections.lasso.z1)))
  
  ####################### COMBINED ##############################
  
  selections.combined = rbind(Stats.z0, Stats.z1) %>% knockoff.filter(fdr=fdr, offset=1) %>% 
    arrange(CHR, Group) %>% 
    arrange(desc(W))  
  print(paste0("Number selections lasso combined: ", nrow(selections.combined)))

  selections.combined = selections.combined %>% left_join(sampled_nonzero_snps, by = c("CHR_Group"))
  
  
  get_false_true_null = function(x) {
    
    x = as.data.frame(t(x))
    
    # can only select a single interaction variable
    if(x$environment.selected == 1) {
      null = all(rowSums(ite[ind.z0, which(map.real.chr21.orig$CHR_Group == x$CHR_Group), drop = FALSE]) == 0)
    } else {
      null = all(rowSums(ite[ind.z1, which(map.real.chr21.orig$CHR_Group == x$CHR_Group), drop = FALSE]) == 0)
    }
    
    return(null)
  }
  
  get_homogeneity = function(x) {
    
    x = as.data.frame(t(x))
    # can only select a single interaction variable
    if(x$environment.selected == 1) {
      homogeneity = mean(rowSums(ite[ind.z0, which(map.real.chr21.orig$CHR_Group == x$CHR_Group), drop = FALSE]) != 0)
      
    } else {
      homogeneity = mean(rowSums(ite[ind.z1, which(map.real.chr21.orig$CHR_Group == x$CHR_Group), drop = FALSE]) != 0)
    }
    
    return(homogeneity)
  }
  
  
  if(nrow(selections.combined) > 0) {
    selections.combined$Null = apply(selections.combined, 1, get_false_true_null)
    selections.combined$homogeneity = apply(selections.combined, 1, get_homogeneity)
    
  } else {
    selections.combined$Null = NA
    selections.combined$homogeneity = NA
    
  }
  
  
  # calculate number of rejections by group / snp
  selections.combined = selections.combined %>% 
    group_by(CHR_Group) %>% 
    mutate(rej_by_group = n()) %>% 
    ungroup()
  
  # calculate fdp
  fdp[b] = mean(selections.combined$Null)
  
  # number of correct snps / groups that are rejected at least once
  unique_correct_rejections = selections.combined %>% 
    filter(!Null) %>% 
    select(CHR_Group) %>% 
    distinct()
  
  sum_double_rejections = selections.combined %>% 
    group_by(CHR_Group) %>% 
    mutate(num_rej =  n()) %>% 
    filter(num_rej > 1) %>% 
    select(CHR_Group) %>% 
    distinct()
  
  power_cond_association[b] = nrow(unique_correct_rejections) /  nrow(sampled_nonzero_snps)
  
  # calculate homogeneity measure: proportion of individuals in the reported subgroup for which the variable 
  # of interest has a non-zero coefficient
  homogeneity[b] = mean(selections.combined$homogeneity)
  
  # calculate power and fdp separately for globally / locally nonzero groups 
  selections.combined = selections.combined %>% 
    mutate(false_rejection_interaction = ifelse((Null & rej_by_group == 1) |
                                                  (rej_by_group == 1 & Interaction == "none") | 
                                                  (rej_by_group == 1 & Interaction != "none" & environment != environment.selected), TRUE, FALSE), 
           rejected_single_environment = ifelse(rej_by_group == 1, TRUE, FALSE))
  
  fdp_interactions[b] = sum(selections.combined$false_rejection_interaction) / sum(selections.combined$rejected_single_environment)
  
  
  power_interactions[b] = nrow(selections.combined %>% 
                                 filter((Interaction != "none" & environment.selected == environment &
                                           rej_by_group == 1) & (Null == FALSE))) / sum(sampled_nonzero_snps$Interaction != "none")
  
  selections.combined = selections.combined %>% 
    mutate(false_rejection_main = ifelse((Null & rej_by_group == 2) | 
                                           (!Null & rej_by_group == 2 & Interaction != "none"), TRUE, FALSE)) 
  
  
  
  fdp_main[b] = length(selections.combined %>% filter(false_rejection_main == TRUE) %>% pull(CHR_Group) %>% unique()) /
    length(selections.combined %>% filter(rej_by_group == 2) %>% pull(CHR_Group) %>% unique())
  
  power_main[b] = length(selections.combined %>% filter((rej_by_group == 2 & Interaction == "none") & (Null == FALSE)) %>% pull(CHR_Group) %>% unique()) / 
    sum(sampled_nonzero_snps$Interaction == "none")
  
  
  power_cond_association_deconstruct_main[b] = length(selections.combined %>% 
                                                        filter((Interaction == "none") & (Null == FALSE)) %>% pull(CHR_Group) %>% unique()) / sum(sampled_nonzero_snps$Interaction == "none")
  
  power_cond_association_deconstruct_int[b] = length(selections.combined %>% 
                                                       filter((Interaction != "none") & (Null == FALSE)) %>% pull(CHR_Group) %>% unique()) / sum(sampled_nonzero_snps$Interaction != "none")
  
  ######## separate metrics z = 0, z = 1 ########
  
  selections.lasso.z0 = selections.lasso.z0 %>% left_join(sampled_nonzero_snps, by = c("CHR_Group"))
  selections.lasso.z0 = selections.lasso.z0 %>% 
    mutate(false_rejection = ifelse((nonzero_all != 1 & nonzero_z10 != 1) | is.na(nonzero), TRUE, FALSE))
  
  
  selections.lasso.z1 = selections.lasso.z1 %>% left_join(sampled_nonzero_snps, by = c("CHR_Group"))
  selections.lasso.z1 = selections.lasso.z1 %>% 
    mutate(false_rejection = ifelse((nonzero_all != 1 & nonzero_z11 != 1) | is.na(nonzero), TRUE, FALSE))
  
  
  # calculate fdp and power
  fdp_z0[b] = mean(selections.lasso.z0$false_rejection)
  power_z0[b] = sum(!is.na(selections.lasso.z0$nonzero)) /  nrow(sampled_nonzero_snps)
  
  power_cond_association_deconstruct_main_z0[b] = length(selections.lasso.z0 %>% 
                                                        filter((Interaction == "none") & (!false_rejection)) %>% pull(CHR_Group) %>% unique()) / sum(sampled_nonzero_snps$Interaction == "none")
  
  power_cond_association_deconstruct_int_z0[b] = length(selections.lasso.z0 %>% 
                                                       filter((Interaction != "none") & (!false_rejection)) %>% pull(CHR_Group) %>% unique()) / sum(sampled_nonzero_snps$Interaction != "none")
  
  homogeneity_z0[b] = 1 - mean(selections.lasso.z0$false_rejection)
  
  
  
  fdp_z1[b] =  mean(selections.lasso.z1$false_rejection)
  power_z1[b] = sum(!is.na(selections.lasso.z1$nonzero)) /  nrow(sampled_nonzero_snps)
  
  power_cond_association_deconstruct_main_z1[b] = length(selections.lasso.z1 %>% 
                                                        filter((Interaction == "none") & (!false_rejection)) %>% pull(CHR_Group) %>% unique()) / sum(sampled_nonzero_snps$Interaction == "none")
  
  power_cond_association_deconstruct_int_z1[b] = length(selections.lasso.z1 %>% 
                                                       filter((Interaction != "none") & (!false_rejection)) %>% pull(CHR_Group) %>% unique()) / sum(sampled_nonzero_snps$Interaction != "none")
  
  homogeneity_z1[b] = 1 - mean(selections.lasso.z1$false_rejection)
  
}

results.df = data.frame(fdp = fdp, 
                        power_cond_association = power_cond_association,
                        homogeneity = homogeneity, 
                        homogeneity_complement = NA, 
                        fdp_interactions = fdp_interactions, 
                        power_interactions = power_interactions, 
                        fdp_main = fdp_main, 
                        power_main = power_main,
                        power_cond_association_deconstruct_main = power_cond_association_deconstruct_main, 
                        power_cond_association_deconstruct_int = power_cond_association_deconstruct_int, 
                        heritability = heritability,
                        sparsity = sparsity, amp = amp, fdr = fdr, propmain = propmain, method = "sepcomb")

print("Column means without NA")
results.df %>% 
  group_by(method) %>%
  summarise_all(~ mean(., na.rm = TRUE))

# Save res
stats.path = paste0("UKB_stats/simint_paper_new_combined_separate_lasso_", 
                    my_population, "_res", resolution, "_sim", "_B", B, "_num_cov", num_covariates)
res.out.file <- sprintf("%s_s%g_a%g_f%g_mp%g_z1%g_z0%g_r%g.txt", stats.path, sparsity, amp, fdr, propmain, prop_z1, prop_z0, random_sample_prop)
results.df %>% write_delim(res.out.file, delim=" ")
cat(sprintf("Results written on: %s\n", res.out.file))



# results for Z = 0
results.df.z0 = data.frame(fdp = fdp_z0, 
                           power_cond_association = power_z0,
                           homogeneity = homogeneity_z0, 
                           homogeneity_complement = NA, 
                           fdp_interactions = NA, 
                           power_interactions = NA, 
                           fdp_main = NA, 
                           power_main = NA,
                           power_cond_association_deconstruct_main = power_cond_association_deconstruct_main_z0, 
                           power_cond_association_deconstruct_int = power_cond_association_deconstruct_int_z0, 
                           heritability = heritability,
                           sparsity = sparsity, amp = amp, fdr = fdr, propmain = propmain, method = "sepz0")

print("Column means without NA")
results.df.z0 %>% 
  group_by(method) %>%
  summarise_all(~ mean(., na.rm = TRUE))

# Save res
stats.path = paste0("UKB_stats/simint_paper_new_z0_separate_lasso_", 
                    my_population, "_res", resolution, "_sim", "_B", B, "_num_cov", num_covariates)
res.out.file <- sprintf("%s_s%g_a%g_f%g_mp%g_z1%g_z0%g_r%g.txt", stats.path, sparsity, amp, fdr, propmain, prop_z1, prop_z0, random_sample_prop)
results.df.z0 %>% write_delim(res.out.file, delim=" ")
cat(sprintf("Results written on: %s\n", res.out.file))

# results for Z = 1
results.df.z1 = data.frame(fdp = fdp_z1, 
                           power_cond_association = power_z1,
                           homogeneity = homogeneity_z1, 
                           homogeneity_complement = NA, 
                           fdp_interactions = NA, 
                           power_interactions = NA, 
                           fdp_main = NA, 
                           power_main = NA,
                           power_cond_association_deconstruct_main = power_cond_association_deconstruct_main_z1, 
                           power_cond_association_deconstruct_int = power_cond_association_deconstruct_int_z1, 
                           heritability = heritability,
                           sparsity = sparsity, amp = amp, fdr = fdr, propmain = propmain, method = "sepz1")

print("Column means without NA")
results.df.z1 %>% 
  group_by(method) %>%
  summarise_all(~ mean(., na.rm = TRUE))

# Save res
stats.path = paste0("UKB_stats/simint_paper_new_z1_separate_lasso_", 
                    my_population, "_res", resolution, "_sim", "_B", B, "_num_cov", num_covariates)
res.out.file <- sprintf("%s_s%g_a%g_f%g_mp%g_z1%g_z0%g_r%g.txt", stats.path, sparsity, amp, fdr, propmain, prop_z1, prop_z0, random_sample_prop)
results.df.z1 %>% write_delim(res.out.file, delim=" ")
cat(sprintf("Results written on: %s\n", res.out.file))




