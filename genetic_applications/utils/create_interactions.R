#!/usr/bin/env Rscript

# UKB GWAS
#
# Class: R script
#
# This script creates interactions between genotypes and covariates

# Load packages
suppressMessages(library(tidyverse))
suppressMessages(library(bigsnpr))

# Input arguments
args <- commandArgs(trailingOnly=TRUE)
my_population <- as.character(args[1])
resolution <- as.numeric(args[2])
covariate_index <- as.integer(args[3])
pheno.name <- as.character(args[4])

print(paste0("Population: ", my_population, 
             " Resolution: ", resolution, 
             " Covariate index: ", covariate_index, 
             " Phenotype: ", pheno.name))


set.seed(2020)
print(pheno.name)
swapped_genotypes_path <- "path_where_swapped_genotypes_are_stored"
print(swapped_genotypes_path)


## Load (screened / nonzero) genotypes ##
covariate_file = paste0(swapped_genotypes_path,"/covariates/covariates_",
                        pheno.name, "_",
                        my_population ,"_sample_all_groups")
covariates = read.csv(covariate_file)

# get last part after "_"
get_last_string <- function(string) {
  parts <- strsplit(string, "_")[[1]]
  return(tail(parts, 1))
}
filtered_col_names <- colnames(covariates)[!grepl("eid", colnames(covariates)) & !grepl("^rev", colnames(covariates))]
covar_names_pretty <- sapply(filtered_col_names, get_last_string)
covar_identifiyer = paste0(sort(covar_names_pretty), collapse = "_")


prescreened.swapped.genotypes.filename <- paste0(swapped_genotypes_path, "chr_merged/prescreened/UKB_swp_", 
                                                 my_population, "_res", resolution, "_" ,pheno.name,
                                                 "_cov_", covar_identifiyer, "_ordered", "_prescreened")
# Attach the "bigSNP" object in R session
cat("Attaching bigSNP object... ")
if(!file.exists(paste0(prescreened.swapped.genotypes.filename, ".bk"))) {
  prescreened.swapped.genotypes.filename.rds = snp_readBed2(paste0(prescreened.swapped.genotypes.filename, ".bed"))
}
obj.bigSNP <- snp_attach(paste0(prescreened.swapped.genotypes.filename, ".rds"))
cat("done.\n")


# Make sure that the rows of the genotypes match the rows of the phenotypes
Subjects <- obj.bigSNP$fam %>% as_tibble()
colnames(Subjects) <- c("FID", "IID", "X1", "X2", "sex", "X3")
Subjects <- Subjects %>% select(FID, IID) %>%
  mutate(FID=as.character(FID), IID=as.character(IID))
if(sum(covariates$eid != Subjects$FID) > 0) {
  stop("Subjects not in correct order!")
}

# remove eid
covariates = covariates %>% dplyr::select(-eid)

#########################
## Create interactions ##
#########################

cat("Creating interactions ... ")
print(dim(obj.bigSNP$genotypes))
print(paste0("NUMBER OF CORES: ", nb_cores()))
G <- obj.bigSNP$genotypes
covariate = covariates[, covariate_index]
head(covariate)
print(length(covariate))
code = obj.bigSNP$genotypes$code256
backingfile.name = paste0(backingfile = "swapped_genotypes_path/tmp/ukb_Gcov_interaction_", 
                          my_population, "_res", resolution)
if(file.exists(paste0(backingfile.name, ".bk"))) {
  file.remove(paste0(backingfile.name, ".bk"))
}
G_cov <- big_copy(G, backingfile = backingfile.name)

big_apply(G_cov, function(G_cov, ind, covariate) {
  # have an idea of progress
  #print(ind[1])
  # access a subset of columns as a standard R matrix
  subset_genotypes = G_cov[,ind]
  # multiply 
  G_cov[, ind] = subset_genotypes*covariate
  # return nothing
  NULL
}, block.size = 500, covariate = covariate, a.combine = 'c')

warnings()
cat("done.\n")

map_interactions = obj.bigSNP$map %>% 
  dplyr::mutate(marker.ID = paste0(marker.ID, ":", colnames(covariates)[covariate_index]))


output <- structure(list(genotypes = G_cov,
                         fam = obj.bigSNP$fam,
                         map = map_interactions),
                    class = "bigSNP")

print("Dimensions:")
print(dim(G_cov))

#########################
## Save interactions ##
#########################

# write as BED/BIM/FAM file triple
interaction_filename = paste0(swapped_genotypes_path,
                              "covariate_interactions/covar_interaction_", 
                              my_population, "_res",
                              resolution, "_cov_", covar_identifiyer, "_cov", covariate_index)

# remove files if they existed
if(covariate_index == 1) {
  to_be_deleted <- list.files(paste0(swapped_genotypes_path, "covariate_interactions/")
                              , pattern = paste0("covar_interaction_", my_population, "_res", resolution))
  
  if(length(to_be_deleted) > 0) {
    file.remove(paste0(swapped_genotypes_path, "covariate_interactions/", to_be_deleted))
  }
  
}

snp_writeBed(output, paste0(interaction_filename, ".bed"))


