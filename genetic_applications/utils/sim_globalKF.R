#!/usr/bin/env Rscript

# Simulation based on UK Biobank genotypes
#
# Class: R script
#
# Simulate performance of vanilla knockoff filter



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

fdp = c()
power_cond_association = c()
homogeneity = c()

power_cond_association_deconstruct_int = c()
power_cond_association_deconstruct_main = c()

heritability = c()


print(paste0("NUMBER ITERATIONS: ", B))

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
  
  # randomly assign each non-zero snps to local and global effects
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
  rm(list=setdiff(ls(), c("y", "covariates", "map.real.chr21.orig.main", "map.real.chr21.main", "map.real.chr21", "map.real.chr21.orig",  "map.real",  "map.real.main", 
                          "G.real.chr21", "G.real.chr21.orig", 
                          "my_population", "sparsity", "fdr", "amp", 
                          "ncores", "dfmax", "B", "b", "resolution",   "num_covariates",
                          "dfmax_lasso_all",  "sample_nonzero_interaction", "sampled_nonzero_snps",
                          "fdp", "power_cond_association", "propmain", "ite", "homogeneity", 
                          "power_main", "fdp_main", "power_interactions", "fdp_interactions",
                          "power_cond_association_deconstruct_int", "power_cond_association_deconstruct_main",
                          "heritability", "prop_z1", "prop_z0", "random_sample_prop")))
  
  # calculate sum of importance statistics from swapped lasso
  W.stats <- function(Z, knockoff) {
    importance <- abs(Z)
    z  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
    zk <- sum(importance[which(knockoff==TRUE)], na.rm=T)
    w <- z-zk
  }
  
  

  ######################### Run lasso + knockoff filter ############################
  
  # Compute scaling factor for the genotypes
  cat("Computing scaling factors for all variants... ")
  scaler <- big_scale()
  G.scale <- scaler(G.real.chr21)
  scaling.factors <- G.scale$scale
  cat("done.\n")
  
  # Run lasso on swapped genotypes for pre-screening 
  cat(sprintf("Fitting sparse linear regression with %d observations, %d variants and NO covariates... ",
              length(y), ncol(G.real.chr21)))
  lasso.fit <- big_spLinReg(G.real.chr21, y.train=y, dfmax=dfmax, ncores=ncores, 
                            covar.train=covariates, pf.covar = rep(0, ncol(covariates)))
  
  # Extract beta from each fold and combine them
  cat("Extracting regression coefficients... ")
  beta <- sapply(1:10, function(k) lasso.fit[[1]][k][[1]]$beta)
  
  # get coefficients
  total_number_genotypes_in_model = nrow(beta) - ncol(covariates)
  beta.variants <- beta[1:total_number_genotypes_in_model,]
  
  # Undo scaling of lasso coefficients
  remaining_columns = attr(lasso.fit, "ind.col") # which genotypes were used for model fitting (index does not include covariates)
  beta.variants <- beta.variants * scaling.factors[remaining_columns]
  
  Beta <- cbind(tibble(CHR=map.real.chr21[remaining_columns,]$CHR,
                       SNP=map.real.chr21[remaining_columns,]$SNP, BP=map.real.chr21[remaining_columns,]$BP, 
                       Knockorr = map.real.chr21[remaining_columns,]$Knockoff),
                as_tibble(beta.variants)) %>% as_tibble()
  colnames(Beta) <- c("CHR", "SNP", "BP", "Knockoff",paste("K", seq(ncol(beta.variants)),sep=""))
  Beta <- Beta %>%
    mutate(Z=(K1+K2+K3+K4+K5+K6+K7+K8+K9+K10)/10,
           Nonzero=(K1!=0)+(K2!=0)+(K3!=0)+(K4!=0)+(K5!=0)+(K6!=0)+(K7!=0)+(K8!=0)+(K9!=0)+(K10!=0)) %>%
    select(CHR, SNP, BP, Knockoff, Z)
  
  Beta <- Beta %>% 
    left_join(map.real.chr21, by = c("CHR", "SNP", "Knockoff", "BP"))
  
  Stats <- Beta %>%
    select("CHR", "Group", "CHR_Group", "SNP", "BP", "Knockoff", "Z") %>%
    group_by(CHR, Group) %>%
    summarize(W = W.stats(abs(Z),Knockoff),
              Lead=which.max(abs(Z)), SNP.lead=SNP[Lead], BP.lead=BP[Lead],
              Size=n()) %>%
    ungroup() %>%
    select(CHR, Group, SNP.lead, BP.lead, Size, W) %>%
    mutate(CHR_Group = paste0(CHR, "_", Group))
  
  # Create the bar plot
  ggplot(Stats, aes(x = reorder(W, - abs(W)), y = W, fill = W > 0)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    scale_fill_manual(values = c("FALSE" = "red", "TRUE" = "blue"), 
                      labels = c("Negative", "Positive")) +
    labs(title = "Stats Vanilla", 
         x = "W (ordered by absolute value)", 
         y = "Value") +
    theme_minimal() +
    theme(axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(),  # Optionally hide x-axis text
          panel.grid.major.x = element_blank(),  # Remove major grid lines for x
          panel.grid.minor.x = element_blank()) 
  
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
  selections.lasso <- Stats %>% knockoff.filter(fdr=fdr, offset=1) %>% 
    arrange(CHR, Group) %>% 
    arrange(desc(W)) 
  
  print(paste0("Number selections lasso: ", nrow(selections.lasso)))
  
  selections.lasso = selections.lasso %>% left_join(sampled_nonzero_snps, by = c("CHR_Group"))
  
  selections.lasso = selections.lasso %>% 
    mutate(false_rejection = ifelse(is.na(nonzero), TRUE, FALSE), 
           correct_rejection = !false_rejection) 
  
  fdp[b] = mean(selections.lasso$false_rejection)
  
  power_cond_association[b] = sum(selections.lasso$correct_rejection) / nrow(sampled_nonzero_snps)
  
  
  get_homogeneity = function(x) {
    
    x = as.data.frame(t(x))
    # vanilla selects all individuals by default always
    homogeneity = mean(rowSums(ite[, which(map.real.chr21.orig$CHR_Group == x$CHR_Group), drop = FALSE]) != 0)
    
    return(homogeneity)
  }
  
  if(nrow(selections.lasso) > 0) {
    selections.lasso$homogeneity = apply(selections.lasso, 1, get_homogeneity)
  } else {
    selections.lasso$homogeneity = NA
  }
  
  
  power_cond_association_deconstruct_main[b] = length(selections.lasso %>% 
                                                        filter((Interaction == "none") & (!false_rejection)) %>% pull(CHR_Group) %>% unique()) / sum(sampled_nonzero_snps$Interaction == "none")
  
  power_cond_association_deconstruct_int[b] = length(selections.lasso %>% 
                                                       filter((Interaction != "none") & (!false_rejection)) %>% pull(CHR_Group) %>% unique()) / sum(sampled_nonzero_snps$Interaction != "none")
  
  
  homogeneity[b] = mean(selections.lasso$homogeneity)
  
}

results.df = data.frame(fdp = fdp, 
                        power_cond_association = power_cond_association,
                        homogeneity = homogeneity, 
                        homogeneity_complement = NA, 
                        fdp_interactions = NA, 
                        power_interactions = NA, 
                        fdp_main = NA, 
                        power_main = NA,
                        power_cond_association_deconstruct_int = power_cond_association_deconstruct_int, 
                        power_cond_association_deconstruct_main = power_cond_association_deconstruct_main, 
                        heritability = heritability,
                        sparsity = sparsity, amp = amp, fdr = fdr, propmain = propmain, method = "sskf")

print("Column means without NA")
results.df %>% 
  group_by(method) %>%
  summarise_all(~ mean(., na.rm = TRUE))

# Save fdp 
stats.path = paste0("UKB_stats/simint_paper_vanilla_lasso_", 
                    my_population, "_res", resolution, "_sim", "_B", B, "_num_cov", num_covariates)
res.out.file <- sprintf("%s_s%g_a%g_f%g_mp%g_z1%g_z0%g_r%g.txt", stats.path, sparsity, amp, fdr, propmain, prop_z1, prop_z0, random_sample_prop)
results.df %>% write_delim(res.out.file, delim=" ")
cat(sprintf("Results written on: %s\n", res.out.file))



