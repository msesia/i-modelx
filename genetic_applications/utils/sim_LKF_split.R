#!/usr/bin/env Rscript


# Simulation based on UK Biobank genotypes
#
# Class: R script
#
# Simulate performance of split analysis (half for train, half for test)

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
num_covariates <- as.numeric(args[12])

library(bigsnpr)
library(tidyverse)

print(version)

my_population = "whitenonbritish"
resolution = 1
num_inter = 1   

fdp = c()
fdp_interactions = c()
power_cond_association = c()
homogeneity = c()
homogeneity_complement = c()

fdp_interactions = c()
power_interactions = c()
fdp_main = c()
power_main = c()

power_cond_association_deconstruct_int = c()
power_cond_association_deconstruct_main = c()

heritability = c()

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

print(paste0("Number of iterations: ", B))

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
  
  # randomly assign each non-zero snp to local or global
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
                          "G.real.chr21", "G.real.chr21.orig",  "my_population", "sparsity", "fdr", "amp", 
                          "ncores", "dfmax", "B", "b", "resolution",  "num_covariates",
                          "num_inter", "dfmax_lasso_all",  "sample_nonzero_interaction", "sampled_nonzero_snps",
                          "fdp", 
                          "fdp_interactions",  "power_cond_association_deconstruct_int", "power_cond_association_deconstruct_main",
                          "power_cond_association", "propmain", "ite", "homogeneity", "homogeneity_complement", 
                          "power_main", "fdp_main", "power_interactions", "fdp_interactions", "heritability", "prop_z1", "prop_z0", "random_sample_prop")))
  
  # calculate sum of importance statistics 
  W.stats <- function(Z, knockoff) {
    importance <- abs(Z)
    z  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
    zk <- sum(importance[which(knockoff==TRUE)], na.rm=T)
    w <- z-zk
  }
  
  W.stats.sum <- function(Z, knockoff) {
    importance <- abs(Z)
    z  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
    zk <- sum(importance[which(knockoff==TRUE)], na.rm=T)
    w <- z+zk
  }

  ######################### load real genotypes and consider half for learning #################################

  index_train = sample(1:nrow(G.real.chr21), round(nrow(G.real.chr21) / 2))
  index_test = seq(1:nrow(G.real.chr21))[-index_train]
  
  G.real.chr21.train = big_copy(G.real.chr21, ind.row = index_train, 
                                backingfile = tempfile(tmpdir = "tmp"))
  
  G.real.chr21.test = big_copy(G.real.chr21, ind.row = index_test, 
                                backingfile = tempfile(tmpdir = "tmp"))
  
  y.train = y[index_train]
  y.test = y[index_test]
  
  covariates.train = covariates[index_train,, drop = FALSE]
  covariates.test = covariates[index_test,, drop = FALSE]

  # Compute scaling factor for the genotypes
  cat("Computing scaling factors for all variants... ")
  scaler <- big_scale()
  G.scale <- scaler(G.real.chr21.train)
  scaling.factors <- G.scale$scale
  cat("done.\n")
  
  # Run lasso on real training genotypes for pre-screening 
  cat(sprintf("Fitting sparse linear regression with %d observations, %d variants and covariates... ",
              length(y), ncol(G.real.chr21.train)))
  lasso.fit <- big_spLinReg(G.real.chr21.train, y.train=y.train, dfmax=dfmax, ncores=ncores, 
                            covar.train=covariates.train, pf.covar = rep(0, ncol(covariates.train)))
  
  # Extract beta from each fold and combine them
  cat("Extracting regression coefficients... ")
  beta <- sapply(1:10, function(k) lasso.fit[[1]][k][[1]]$beta)
  
  # extract coefficients
  total_number_genotypes_in_model = nrow(beta) - ncol(covariates.train)
  beta.variants <- beta[1:total_number_genotypes_in_model,]
  
  # Undo scaling of lasso coefficients
  remaining_columns = attr(lasso.fit, "ind.col") # which genotypes were used for model fitting (index does not include covariates)
  beta.variants <- beta.variants * scaling.factors[remaining_columns]
  Beta <- cbind(tibble(CHR=map.real.chr21[remaining_columns, ]$CHR,
                       SNP=map.real.chr21[remaining_columns, ]$SNP,
                       BP=map.real.chr21[remaining_columns, ]$BP, 
                       Knockoff = map.real.chr21[remaining_columns, ]$Knockoff),
                as_tibble(beta.variants)) %>% as_tibble()
  colnames(Beta) <- c("CHR", "SNP", "BP", "Knockoff",paste("K", seq(ncol(beta.variants)),sep=""))
  Beta <- Beta %>%
    mutate(Z=(K1+K2+K3+K4+K5+K6+K7+K8+K9+K10)/10,
           Nonzero=(K1!=0)+(K2!=0)+(K3!=0)+(K4!=0)+(K5!=0)+(K6!=0)+(K7!=0)+(K8!=0)+(K9!=0)+(K10!=0)) %>%
    select(CHR, SNP, BP, Knockoff, Z)
  
  Beta <- Beta %>% 
    left_join(map.real.chr21, by = c("CHR", "SNP", "BP","Knockoff"))
  
  Stats <- Beta %>%
    select("CHR", "Group", "CHR_Group", "SNP", "BP", "Knockoff", "Z") %>%
    group_by(CHR, Group) %>%
    summarize(W.sum = W.stats(abs(Z),Knockoff),
              Lead=which.max(abs(Z)), SNP.lead=SNP[Lead], BP.lead=BP[Lead],
              Size=n()) %>%
    ungroup() %>%
    select(CHR, Group, SNP.lead, BP.lead, Size, W.sum) %>%
    mutate(CHR_Group = paste0(CHR, "_", Group))
  
  # select nonzero snps
  nonzero.screened.snp =  Beta %>%
    select("CHR", "Group", "CHR_Group", "SNP", "BP", "Knockoff", "Z") %>%
    group_by(SNP) %>%
    summarize(W = W.stats(abs(Z),Knockoff)) %>%
    ungroup() %>%
    select(SNP, W) %>% 
    filter(W > 0) %>% 
    mutate(nonzero = 1)
  
  # save nonzero snps 
  map.real.chr21 = map.real.chr21 %>% 
    left_join(nonzero.screened.snp, by = c("SNP"))
  
  index.nonzero.snps = which(map.real.chr21$nonzero == 1)
  
  map.real.chr21.nonzero.snps = map.real.chr21[index.nonzero.snps, ] 
  
  
  ################## lasso with interactions - still on train #######
  
  # get swapped genotypes for nonzero
  G.real.chr21.train.nonzero = as_FBM(G.real.chr21.train[, index.nonzero.snps],
                                   backingfile = tempfile(tmpdir = "tmp"))
  
  # create interactions
  
  if(num_covariates == 1) {
    G.real.chr21.train.nonzero.int = matrix(NA, nrow = nrow(G.real.chr21.train.nonzero), ncol = ncol(G.real.chr21.train.nonzero)*3)
  }
  
  if(num_covariates == 2) {
    G.real.chr21.train.nonzero.int = matrix(NA, nrow = nrow(G.real.chr21.train.nonzero), ncol = ncol(G.real.chr21.train.nonzero)*5)
  }
  
  G.real.chr21.train.nonzero.int[, 1:ncol(G.real.chr21.train.nonzero)] = G.real.chr21.train.nonzero[]
  
  G.real.chr21.train.nonzero.int[, (ncol(G.real.chr21.train.nonzero)+1):(2*ncol(G.real.chr21.train.nonzero))] = G.real.chr21.train.nonzero[]*covariates.train[, 1]
  G.real.chr21.train.nonzero.int[, (ncol(G.real.chr21.train.nonzero)*2+1):(3*ncol(G.real.chr21.train.nonzero))] = G.real.chr21.train.nonzero[]*(1-covariates.train[, 1])
  
  
  if(num_covariates == 2) {
    G.real.chr21.train.nonzero.int[, (ncol(G.real.chr21.train.nonzero)*3+1):(4*ncol(G.real.chr21.train.nonzero))] = G.real.chr21.train.nonzero[]*covariates.train[, 2]
    G.real.chr21.train.nonzero.int[, (ncol(G.real.chr21.train.nonzero)*4+1):(5*ncol(G.real.chr21.train.nonzero))] = G.real.chr21.train.nonzero[]*(1-covariates.train[, 2])
  }
  
  
  G.real.chr21.train.nonzero.int = as_FBM(G.real.chr21.train.nonzero.int,
                                       backingfile = tempfile(tmpdir = "tmp"))
  
  # create map for interactions 
  map.real.chr21.nonzero.snps.main = map.real.chr21.nonzero.snps %>% mutate(Interaction = "none")
  map.real.chr21.nonzero.snps.intz1 = map.real.chr21.nonzero.snps %>% mutate(Interaction = "Z1")
  map.real.chr21.nonzero.snps.intz1rev = map.real.chr21.nonzero.snps %>% mutate(Interaction = "rev_Z1")
  
  if(num_covariates == 1) {
    map.real.chr21.nonzero.snps.inter = rbind(map.real.chr21.nonzero.snps.main, 
                                               map.real.chr21.nonzero.snps.intz1, 
                                               map.real.chr21.nonzero.snps.intz1rev) 
  }
  
  if(num_covariates == 2) {
    map.real.chr21.nonzero.snps.intz2 = map.real.chr21.nonzero.snps %>% mutate(Interaction = "Z2")
    map.real.chr21.nonzero.snps.intz2rev = map.real.chr21.nonzero.snps %>% mutate(Interaction = "rev_Z2")
    map.real.chr21.nonzero.snps.inter = rbind(map.real.chr21.nonzero.snps.main, 
                                               map.real.chr21.nonzero.snps.intz1, 
                                               map.real.chr21.nonzero.snps.intz1rev, 
                                               map.real.chr21.nonzero.snps.intz2, 
                                               map.real.chr21.nonzero.snps.intz2rev) 
  }
  
  
  # run lasso on data with interaction
  
  # Compute scaling factor for the genotypes
  cat("Computing scaling factors for all variants... ")
  scaler.inter <- big_scale()
  G.scale.inter <- scaler.inter(G.real.chr21.train.nonzero.int)
  scaling.factors.inter <- G.scale.inter$scale
  cat("done.\n")
  
  # Run lasso on swapped genotypes for pre-screening 
  cat(sprintf("Fitting sparse linear regression with %d observations, %d variants and covariates... ",
              length(y), ncol(G.real.chr21.train.nonzero.int)))
  lasso.fit.inter <- big_spLinReg(G.real.chr21.train.nonzero.int, y.train=y.train, dfmax=dfmax, ncores=ncores, 
                                  covar.train=covariates.train, pf.covar = rep(0, ncol(covariates.train)))
  
  # Extract beta from each fold and combine them
  cat("Extracting regression coefficients... ")
  beta.inter <- sapply(1:10, function(k) lasso.fit.inter[[1]][k][[1]]$beta)
  
  # extract coefficients 
  total_number_genotypes_in_model = nrow(beta.inter) - ncol(covariates)
  beta.variants.inter <- beta.inter[1:total_number_genotypes_in_model,]
  
  # Undo scaling of lasso coefficients
  remaining_columns = attr(lasso.fit.inter, "ind.col") # which genotypes were used for model fitting (index does not include covariates)
  beta.variants.inter <- beta.variants.inter * scaling.factors.inter[remaining_columns]
  Beta.inter <- cbind(tibble(CHR=map.real.chr21.nonzero.snps.inter[remaining_columns,]$CHR,
                             SNP=map.real.chr21.nonzero.snps.inter[remaining_columns,]$SNP, 
                             BP=map.real.chr21.nonzero.snps.inter[remaining_columns,]$BP, 
                             Knockoff = map.real.chr21.nonzero.snps.inter[remaining_columns,]$Knockoff, 
                             Interaction = map.real.chr21.nonzero.snps.inter[remaining_columns,]$Interaction),
                      as_tibble(beta.variants.inter)) %>% as_tibble()
  colnames(Beta.inter) <- c("CHR", "SNP", "BP", "Knockoff","Interaction",paste("K", seq(ncol(beta.variants)),sep=""))
  Beta.inter <- Beta.inter %>%
    mutate(Z=(K1+K2+K3+K4+K5+K6+K7+K8+K9+K10)/10,
           Nonzero=(K1!=0)+(K2!=0)+(K3!=0)+(K4!=0)+(K5!=0)+(K6!=0)+(K7!=0)+(K8!=0)+(K9!=0)+(K10!=0)) %>%
    select(CHR, SNP, BP, Knockoff, Interaction,Z)
  
  
  # Extract the estimated coefficients
  Lasso.res.inter <- Beta.inter %>%
    inner_join(map.real.chr21.nonzero.snps.inter, by = c("CHR", "SNP", "BP", "Knockoff", "Interaction")) %>%
    select(CHR, SNP, BP, Z, Group, Knockoff, Interaction)
  
  # choose top interactions using knockoff importance statistics (difference! of abs values, not sum, not cloaked data)
  
  Stats.inter <- Lasso.res.inter %>%
    select("CHR", "Group", "SNP", "BP", "Knockoff", "Z", "Interaction") %>%
    group_by(CHR, Group, Interaction) %>%
    summarize(W = W.stats(abs(Z),Knockoff),
              Lead=which.max(abs(Z)), SNP.lead=SNP[Lead], BP.lead=BP[Lead],
              Size=n()) %>%
    ungroup() %>%
    select(CHR, Group, SNP.lead, BP.lead, Size, W, Interaction)
  
  Stats.inter = Stats.inter %>% mutate(CHR_Group = paste0(CHR, "_", Group)) %>% 
    rename(Interaction_orig = Interaction) %>%
    mutate(Interaction = sub("^rev_", "", Interaction_orig))

  # choose top W by group
  partition_tmp <- Stats.inter %>%
    filter(W > 0) %>% # then we can filter for positive ones!
    group_by(CHR, Group) %>% #maybe change later: choose by chr, group, SNP.raw
    dplyr::mutate(number_unique_interactions = n_distinct(Interaction))  %>% 
    arrange(CHR, Group) %>% 
    dplyr::filter(!(number_unique_interactions > 1 & Interaction == "none")) %>% 
    # if there are two interactions picked (e.g. with sex and rev_sex) pick only the top one
    group_by(CHR, Group, Interaction) %>% 
    top_n(pmin(n(),1), W) 
  
  partition = partition_tmp %>% 
    group_by(CHR, Group) %>%
    top_n(pmin(n(),num_inter), W) %>% 
    dplyr::select(-number_unique_interactions)
  
  chr_group_with_main_effect = Stats.inter %>% filter(W > 0) %>%
    filter(Interaction == "none") %>% 
    dplyr::select(CHR, Group) %>% 
    unique() %>% 
    mutate(has_main_effect = 1)
  
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
  partition = partition %>% mutate(CHR_Group = paste0(CHR, "_", Group))
  unique_chr_group = unique(partition$CHR_Group)
  all_interaction_vars = list()
  for(r in unique_chr_group) {
    # filter partition to the group; arrange ensures that the ordering for all is the same
    chr_group_interaction = partition %>% dplyr::filter(CHR_Group == r) %>% arrange(Interaction)
    all_interaction_vars = c(all_interaction_vars, list(chr_group_interaction$Interaction))
  }
  
  
  ################ lasso env specific - now on test  ################
  
  unique_interaction_vars = unique(all_interaction_vars)
  all_environments = lapply(unique_interaction_vars, function(x) partition_covariates(data.frame(covariates.test), x))
  
 
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
  
  all_stats_lasso_all = c()
  
  for(i in 1:length(unique_interaction_vars)) {
    
    print(paste0("Working on interaction ", unique_interaction_vars[[i]]))
    
    # get all variables with this interaction combination
    relevant_chr_group = chr_group_each_interaction[[i]]
  
    # get environments for this interaction
    environments = all_environments[[i]]
    unique_environments = sort(unique(environments))
    
    
    for(e in unique_environments) {
      
      print(paste0("Working on environment ", e, " out of ", length(unique_environments)))
      
      index_e = which(environments == e)
      index_chr_group = which(map.real.chr21$CHR_Group %in% relevant_chr_group)
      
      G.real.subset.env = big_copy(G.real.chr21.test, ind.row = index_e,
                                   backingfile = tempfile(tmpdir = "tmp"))
      
      
      # run lasso on this subset
      scaler <- big_scale()
      G.scale <- scaler(G.real.subset.env)
      scaling.factors <- G.scale$scale
      
      # Run lasso 
      covariates_environment = covariates.test[index_e, , drop = FALSE]
      
      # filter out covariates that do not have any variance
      variance_each_covar = sapply(as.data.frame(covariates_environment), function(x) c(var=var(x)))
      print("null variance covars ")
      variance_each_covar[which(variance_each_covar == 0)]
      names(variance_each_covar) = NULL
      
      if(sum(variance_each_covar == 0) > 0) {
        covariates_environment = as.matrix(covariates_environment[, -which(variance_each_covar == 0)])
      } else {
        covariates_environment = as.matrix(covariates_environment)
      }
      
      
      # Run lasso 
      lasso.fit <- big_spLinReg(G.real.subset.env, y.train=y.test[index_e], dfmax=dfmax_lasso_all, ncores=ncores, 
                                covar.train=covariates_environment, pf.covar = rep(0, ncol(covariates_environment))) #, 
      
      # Extract beta from each fold and combine them
      beta <- sapply(1:10, function(k) lasso.fit[[1]][k][[1]]$beta)
      
      # get coefficients (no covariates)
      total_number_genotypes_in_model = nrow(beta) - ncol(covariates_environment)
      beta.variants <- beta[1:total_number_genotypes_in_model,]
      
      # Undo scaling of lasso coefficients
      remaining_columns = attr(lasso.fit, "ind.col") # which genotypes were used for model fitting
      beta.variants <- beta.variants * scaling.factors[remaining_columns]
      Beta <- cbind(tibble(CHR=map.real.chr21[remaining_columns, ]$CHR,
                           SNP=map.real.chr21[remaining_columns, ]$SNP,
                           BP=map.real.chr21[remaining_columns, ]$BP, 
                           Knockoff=map.real.chr21[remaining_columns, ]$Knockoff),
                    as_tibble(beta.variants)) %>% as_tibble()
      colnames(Beta) <- c("CHR", "SNP", "BP", "Knockoff",paste("K", seq(ncol(beta.variants)),sep=""))
      Beta <- Beta %>%
        mutate(Z=(K1+K2+K3+K4+K5+K6+K7+K8+K9+K10)/10,
               Nonzero=(K1!=0)+(K2!=0)+(K3!=0)+(K4!=0)+(K5!=0)+(K6!=0)+(K7!=0)+(K8!=0)+(K9!=0)+(K10!=0)) %>%
        select(CHR, SNP, BP, Knockoff,Z)
      
      # Extract the estimated coefficients
      Stats <- Beta %>%
        inner_join(map.real.chr21, by = c("CHR", "SNP", "BP", "Knockoff")) %>%
        select("CHR", "Group", "SNP", "BP", "Knockoff", "Z") %>%
        group_by(CHR, Group) %>%
        summarize(W = W.stats(abs(Z),Knockoff),
                  Lead=which.max(abs(Z)), SNP.lead=SNP[Lead], BP.lead=BP[Lead],
                  Size=n()) %>%
        ungroup() %>%
        select(CHR, Group, SNP.lead, BP.lead, Size, W) %>% 
        mutate(Interaction = paste0(unique_interaction_vars[[i]], collapse = ", "), 
               environment = e)
      
      Stats_filtered = Stats %>% 
        mutate(CHR_Group = paste0(CHR, "_", Group)) %>% 
        filter(CHR_Group %in% relevant_chr_group)
      
      all_stats_lasso_all = rbind(all_stats_lasso_all, Stats_filtered)
      
    }
    
    
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
  
  
  # lasso results with interaction
  selections.lasso <- all_stats_lasso_all %>% knockoff.filter(fdr=fdr, offset=1) %>% 
    arrange(CHR, Group) %>% 
    arrange(desc(W)) %>% 
    mutate(CHR_Group = paste0(CHR, "_", Group))
  
  selections.lasso = selections.lasso %>% rename(Interaction.selected = Interaction, 
                                                 environment.selected = environment)
  
  print(paste0("Number selections lasso: ", nrow(selections.lasso)))
  
  selections.lasso = selections.lasso %>% 
    left_join(sampled_nonzero_snps, by = c("CHR_Group")) %>% 
    relocate(Interaction, .after = Interaction.selected) %>% 
    relocate(environment, .after = environment.selected) %>% 
    relocate(nonzero, .after = W)
  
  get_false_true_null = function(x) {
    
    x = as.data.frame(t(x))
    interactions_selected = sort(unlist(strsplit(x$Interaction.selected, ", ")))
    
    index_interaction = which(sapply(unique_interaction_vars, function(x) identical(x, interactions_selected)))
    
    index_interaction_environment = which(all_environments[[index_interaction]] == x$environment.selected)
    
    null = all(rowSums(ite[index_interaction_environment, which(map.real.chr21.orig$CHR_Group == x$CHR_Group), drop = FALSE]) == 0)
    
    return(null)
  }
  
  get_homogeneity = function(x) {
    
    x = as.data.frame(t(x))
    interactions_selected = sort(unlist(strsplit(x$Interaction.selected, ", ")))
    
    index_interaction = which(sapply(unique_interaction_vars, function(x) identical(x, interactions_selected)))
    
    index_interaction_environment = which(all_environments[[index_interaction]] == x$environment.selected)
    
    homogeneity = mean(rowSums(ite[index_interaction_environment, which(map.real.chr21.orig$CHR_Group == x$CHR_Group), drop = FALSE]) != 0)
    
    return(homogeneity)
  }
  
  
  get_homogeneity_complement = function(x) {
    
    x = as.data.frame(t(x))
    interactions_selected = sort(unlist(strsplit(x$Interaction.selected, ", ")))
    
    index_interaction = which(sapply(unique_interaction_vars, function(x) identical(x, interactions_selected)))
    
    index_not_interaction_environment = which(all_environments[[index_interaction]] != x$environment.selected)
    
    homogeneity_complement = mean(rowSums(ite[index_not_interaction_environment, which(map.real.chr21.orig$CHR_Group == x$CHR_Group), drop = FALSE]) != 0)
    
    return(homogeneity_complement)
  }
  
  
  if(nrow(selections.lasso) > 0) {
    selections.lasso$Null = apply(selections.lasso, 1, get_false_true_null)
    selections.lasso$homogeneity = apply(selections.lasso, 1, get_homogeneity)
    
    selections.lasso$homogeneity_complement = apply(selections.lasso, 1, get_homogeneity_complement)
  } else {
    selections.lasso$Null = NA
    selections.lasso$homogeneity = NA
    selections.lasso$homogeneity_complement = NA
    
  }
  
  # calculate number of rejections by group / snp
  selections.lasso = selections.lasso %>% 
    group_by(CHR_Group) %>% 
    mutate(rej_by_group = n()) %>% 
    ungroup()
  
  # calculate fdp
  fdp[b] = mean(selections.lasso$Null)
  
  # number of correct snps / groups that are rejected at least once
  unique_correct_rejections = selections.lasso %>% 
    filter(!Null) %>% 
    select(CHR_Group) %>% 
    distinct()
  
  power_cond_association[b] = nrow(unique_correct_rejections) /  nrow(sampled_nonzero_snps)
  
  # calculate homogeneity measure: proportion of individuals in the reported subgroup for which the variable 
  # of interest has a non-zero coefficient
  homogeneity[b] = mean(selections.lasso$homogeneity)
  
  homogeneity_complement[b] = mean(selections.lasso$homogeneity_complement, na.rm = TRUE)
  
  # calculate power and fdp separately for globally / locally nonzero groups 
  
  # calculate "false discovery" for interactions
  selections.lasso = selections.lasso %>% 
    mutate(false_rejection_interaction = ifelse((Null & Interaction.selected != "none" & rej_by_group == 1 ) | 
                                                  (!Null & Interaction == "none" & Interaction.selected != "none" & rej_by_group == 1) |
                                                  (!Null & Interaction != "none" & Interaction.selected != "none" & rej_by_group == 1 &
                                                     environment.selected != environment),
                                                TRUE, FALSE), 
           rejected_single_environment = ifelse(Interaction.selected != "none" & rej_by_group == 1, TRUE, FALSE))
  
  fdp_interactions[b] = sum(selections.lasso$false_rejection_interaction) / sum(selections.lasso$rejected_single_environment)
  
  
  power_interactions[b] = nrow(selections.lasso %>% 
                                 filter((Interaction != "none" & 
                                           Interaction.selected == Interaction & 
                                           environment.selected == environment &
                                           rej_by_group == 1) & (Null == FALSE))) / sum(sampled_nonzero_snps$Interaction != "none")
  
  selections.lasso = selections.lasso %>% 
    mutate(false_rejection_main = ifelse((Null & Interaction.selected == "none") | 
                                           (!Null & Interaction.selected == "none" & Interaction != "none"), TRUE, FALSE)) 
  
  
  
  fdp_main[b] = sum(selections.lasso$false_rejection_main) / nrow(selections.lasso %>% filter(Interaction.selected == "none"))
  power_main[b] = nrow(selections.lasso %>% filter((Interaction.selected == "none" & Interaction == "none") & (Null == FALSE))) / sum(sampled_nonzero_snps$Interaction == "none")
  
  power_cond_association_deconstruct_main[b] = length(selections.lasso %>% 
                                                        filter((Interaction == "none") & (Null == FALSE)) %>% pull(CHR_Group) %>% unique()) / sum(sampled_nonzero_snps$Interaction == "none")
  
  power_cond_association_deconstruct_int[b] = length(selections.lasso %>% 
                                                       filter((Interaction != "none") & (Null == FALSE)) %>% pull(CHR_Group) %>% unique()) / sum(sampled_nonzero_snps$Interaction != "none")
  
  
  
}

results.df = data.frame(fdp = fdp, 
                        power_cond_association = power_cond_association,
                        homogeneity = homogeneity, 
                        homogeneity_complement = homogeneity_complement, 
                        fdp_interactions = fdp_interactions, 
                        power_interactions = power_interactions, 
                        fdp_main = fdp_main, 
                        power_main = power_main,
                        power_cond_association_deconstruct_main = power_cond_association_deconstruct_main, 
                        power_cond_association_deconstruct_int = power_cond_association_deconstruct_int, 
                        heritability = heritability,
                        sparsity = sparsity, amp = amp, fdr = fdr, propmain = propmain, method = "split")

print("Column means without NA")
results.df %>% 
  group_by(method) %>%
  summarise_all(~ mean(., na.rm = TRUE))

# Save fdp 
stats.path = paste0("UKB_stats/simint_paper_split_lasso_",
                    my_population, "_res", resolution, "_sim", "_B", B, "_num_cov", num_covariates)
res.out.file <- sprintf("%s_s%g_a%g_f%g_mp%g_z1%g_z0%g_r%g.txt", stats.path, sparsity, amp, fdr, propmain, prop_z1, prop_z0, random_sample_prop)
results.df %>% write_delim(res.out.file, delim=" ")
cat(sprintf("Results written on: %s\n", res.out.file))



