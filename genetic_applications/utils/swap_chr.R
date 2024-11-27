#!/usr/bin/env Rscript

# UKB Gwas 
#
# Class: R script
#
# Create cloaked data for a specific chromosome (randomly shuffled knockoff and real genotypes)

# Load packages
suppressMessages(library(tidyverse))
suppressMessages(library(bigsnpr))

args <- commandArgs(trailingOnly=TRUE)
my_population <- as.character(args[1])
resolution  <- as.integer(args[2])
ncores  <- as.integer(args[3])
chromosome  <- as.integer(args[4])

# set paths
outputpath = paste0("swapped_genotypes/chr/")
temp_path = "swapped_genotypes/tmp/"
group_path = "group_information/"

# get WHR files
bed_bim_fam_filename <- paste0("real_data_file")

set.seed(2024)

######## get the real genotypes #########

cat("Reading genotypes ... ")
print(bed_bim_fam_filename)

plinkfile = paste0(bed_bim_fam_filename,
                   ".bed")
if(!file.exists(paste0(bed_bim_fam_filename, ".bk"))) {
  x = snp_readBed2(plinkfile)
}
cat("done.\n")


# Attach the "bigSNP" object in R session
cat("Attaching bigSNP object... ")
obj.bigSNP.real.full <- snp_attach(paste0(bed_bim_fam_filename,".rds"))
cat("done.\n")


backingfilename_orig_subset_chr = paste0(temp_path, "chr_spec_genotype_real_", 
                                         "_", my_population, "_res", resolution, "_chr", chromosome)

if(file.exists(paste0(backingfilename_orig_subset_chr, ".bk"))) {
  file.remove(paste0(backingfilename_orig_subset_chr, ".bk"))
}

if(file.exists(paste0(backingfilename_orig_subset_chr, ".rds"))) {
  file.remove(paste0(backingfilename_orig_subset_chr, ".rds"))
}

obj.bigSNP.real.chr.file <- snp_subset(obj.bigSNP.real.full,
                                       ind.col = which(obj.bigSNP.real.full$map$chromosome == chromosome),
                                       backingfile = paste0(backingfilename_orig_subset_chr))
obj.bigSNP.real <- snp_attach(paste0(obj.bigSNP.real.chr.file))

# Extract list of variants
map.real.orig <- obj.bigSNP.real$map %>% as_tibble()
map.real <- map.real.orig
colnames(map.real) <- c("CHR", "SNP", "gd", "BP", "a1", "a2")
map.real <- map.real %>% dplyr::select(CHR, SNP, BP)

# Get aliases for useful slots
G.real   <- obj.bigSNP.real$genotypes

map.real = map.real %>% mutate(Knockoff = ifelse(endsWith(SNP, ".k"), TRUE, FALSE), 
                               SNP = sub(".k$", "", SNP)) %>% 
  dplyr::select(CHR, SNP, BP, Knockoff)

# Import Group Information 
grp.file <- c()
chr <- seq(1, 22)
for(c in chr) {
  
  myfilename <- paste0(group_path, "ukb_gen_chr", c, "_ibd1_res", resolution, "_grp",".txt")
  myfile <- read.csv(myfilename, sep="") 
  
  myfile$CHR <- c
  
  grp.file <- rbind(grp.file, myfile)
  
}
# from https://rdrr.io/cran/matchingR/src/R/utils.R
repcol <- function(x, n) {
  s <- NCOL(x)
  if (length(n) == 1) {
    return(matrix(x[, rep(1:s, each = n)], nrow = NROW(x), ncol = s * n))
  }
  matrix(x[, rep(1:s, n)], nrow = NROW(x), ncol = sum(n))
}

map.real <- map.real %>% 
  mutate(CHR = as.integer(CHR)) %>%
  left_join(grp.file, by = c("SNP", "CHR")) %>% 
  mutate(index_for_apply = Group + 1) %>% 
  mutate(CHR_Group =paste0(CHR, "_", Group)) %>% 
  mutate(CHR_Index =paste0(CHR, "_", index_for_apply))

ind.real.orig = which(map.real$Knockoff == FALSE)
ind.real.knock = which(map.real$Knockoff == TRUE)

backingfilename_orig = paste0(temp_path, "subset_genotypes_orig", 
                              "_", my_population, "_res", resolution, "_chr", chromosome)

cat("subset bigSNP object... ")
if(file.exists(paste0(backingfilename_orig, ".bk"))) {
  file.remove(paste0(backingfilename_orig, ".bk"))
}

if(file.exists(paste0(backingfilename_orig, ".rds"))) {
  file.remove(paste0(backingfilename_orig, ".rds"))
}

obj.bigSNP.real.orig.file <- snp_subset(obj.bigSNP.real,
                                        ind.col = ind.real.orig,
                                        backingfile = paste0(backingfilename_orig))
obj.bigSNP.real.orig <- snp_attach(paste0(obj.bigSNP.real.orig.file))
G.real.original <- obj.bigSNP.real.orig$genotypes
cat("done.\n")

# same for knockoffs 
backingfilename_knock = paste0(temp_path, "subset_genotypes_knock", 
                               "_", my_population, "_res", resolution, "_chr", chromosome)

cat("subset bigSNP object for knockoffs ... ")
if(file.exists(paste0(backingfilename_knock, ".bk"))) {
  file.remove(paste0(backingfilename_knock, ".bk"))
}

if(file.exists(paste0(backingfilename_knock, ".rds"))) {
  file.remove(paste0(backingfilename_knock, ".rds"))
}

obj.bigSNP.real.knock.file <- snp_subset(obj.bigSNP.real,
                                         ind.col = ind.real.knock,
                                         backingfile = paste0(backingfilename_knock))
obj.bigSNP.real.knock <- snp_attach(paste0(obj.bigSNP.real.knock.file))
G.real.knock <- obj.bigSNP.real.knock$genotypes
cat("done.\n")

map.real.original = map.real[ind.real.orig, ]
map.real.knock = map.real[ind.real.knock, ]

index_apply = unique(map.real.original$CHR_Group)

n = nrow(G.real)

# create matrices holding genotypes for swapping
backingfilename_swapped_orig = paste0(temp_path, "subset_genotypes_orig_swapped", 
                                      "_", my_population, "_res", resolution, "_chr", chromosome)

if(file.exists(paste0(backingfilename_swapped_orig, ".bk"))) {
  file.remove(paste0(backingfilename_swapped_orig, ".bk"))
}

if(file.exists(paste0(backingfilename_swapped_orig, ".rds"))) {
  file.remove(paste0(backingfilename_swapped_orig, ".rds"))
}

G.real.original.swapped <- FBM(nrow(G.real.original), ncol(G.real.original),
                               type = "raw", 
                               backingfile = backingfilename_swapped_orig)

G.real.original.swapped <- add_code256(G.real.original.swapped, code = G.real.original$code256)



backingfilename_swapped_knock = paste0(temp_path, "subset_genotypes_knock_swapped", 
                                       "_", my_population, "_res", resolution, "_chr", chromosome)

if(file.exists(paste0(backingfilename_swapped_knock, ".bk"))) {
  file.remove(paste0(backingfilename_swapped_knock, ".bk"))
}

if(file.exists(paste0(backingfilename_swapped_knock, ".rds"))) {
  file.remove(paste0(backingfilename_swapped_knock, ".rds"))
}

G.real.knock.swapped <- FBM(nrow(G.real.knock), ncol(G.real.knock),
                            type = "raw", 
                            backingfile = backingfilename_swapped_knock)

G.real.knock.swapped <- add_code256(G.real.knock.swapped, code = G.real.knock$code256)


triple_filename = paste0(outputpath, 
                         "UKB_swp_", 
                         my_population, 
                         "_chr", 
                         chromosome, 
                         "_res", 
                         resolution, 
                         "_joint_whr")

if(file.exists(paste0(triple_filename, "_vswap", ".bk"))) {
  file.remove(paste0(triple_filename, "_vswap", ".bk"))
}

# create empty matrix for swapping indices
V.swap.matrix <- FBM(nrow(G.real.original.swapped), 
                     ncol(G.real.original.swapped), 
                     backingfile = paste0(triple_filename, "_vswap"))

set.seed(2024)

#### SWAP  ##
big_apply(G.real.original.swapped, function(G.real.original.swapped,
                                            G.real.knock.swapped, ind, n, V.swap.matrix, map.real.original, map.real.knock) {

  # get index of chr_group
  index.genotypes.original = which(map.real.original$CHR_Group == ind[1])
  index.genotypes.knock = which(map.real.knock$CHR_Group == ind[1])
  
  # subset genotypes
  G.subset.real.original <- G.real.original[,index.genotypes.original, drop = FALSE]
  G.subset.real.knock <- G.real.knock[,index.genotypes.knock, drop = FALSE]
  
  # randomly sample swapping matrix
  V.swap.single.draw = matrix(rbinom(n*1,1,1/2), n)
  V.swap = repcol(V.swap.single.draw, ncol(G.subset.real.original))
  
  V.swap.matrix[,index.genotypes.original] <- V.swap
  
  G.real.original.swapped[,index.genotypes.original] = G.subset.real.original*(1-V.swap)+G.subset.real.knock*(V.swap)
  G.real.knock.swapped[,index.genotypes.knock] = G.subset.real.original*(V.swap)+G.subset.real.knock*(1-V.swap)
  
  # return nothing
  NULL
  
}, block.size = 1, a.combine = 'cbind',
ind = index_apply,
G.real.knock.swapped = G.real.knock.swapped, map.real.knock = map.real.knock,
n = n, ncores = 1, V.swap.matrix = V.swap.matrix, map.real.original = map.real.original)


# compare to manual with v.swap.matrix - only works in resolution 0!
#X.out.1 <- G.real.original[]*(1-V.swap.matrix[])+G.real.knock[]*(V.swap.matrix[])
#X.out.2 <- G.real.original[]*(V.swap.matrix[])+G.real.knock[]*(1-V.swap.matrix[])
#X.out <- cbind(X.out.1, X.out.2)

#all.equal(X.out.1, G.real.original.swapped[])
#all.equal(X.out.2, G.real.knock.swapped[])

### COMBINE ####

if(file.exists(paste0(triple_filename, "_combined", ".bk"))) {
  file.remove(paste0(triple_filename, "_combined", ".bk"))
}

G.combined = FBM(nrow(G.real.original.swapped), 2*ncol(G.real.original.swapped), type = "raw", 
                 backingfile = paste0(triple_filename, "_combined"))

G.combined = add_code256(G.combined, G.real.original$code256)

big_apply(G.real.original.swapped, function(G.real.original.swapped, G.combined, G.real.knock.swapped, p, ind) {
  # have an idea of progress
  
  if(ind[1] %% 1000 == 0) {
    print(ind[1])
  }
  
  G.combined[, ind] <- G.real.original.swapped[, ind]
  G.combined[, ind + p] <- G.real.knock.swapped[, ind]
  
  # return nothing
  NULL
  
}, block.size = 1, a.combine = 'cbind',
ind = cols_along(G.real.original.swapped),
G.real.knock.swapped = G.real.knock.swapped,
G.combined = G.combined,
p = ncol(G.real.original.swapped),
ncores = 1)

# check results
#all.equal(G.combined[, 1:ncol(G.real.original.swapped)][], G.real.original.swapped[])
#all.equal(G.combined[, (ncol(G.real.original.swapped) + 1):(2*ncol(G.real.original.swapped))][], G.real.knock.swapped[])

#### REORDER ##### 
map.swapped = rbind(map.real.original %>% mutate(SNP_swap = paste0(SNP, ".a"), 
                                                 SNP_knock = SNP),
                    map.real.knock %>% mutate(SNP_swap = paste0(SNP, ".b"), 
                                              SNP_knock = paste0(SNP, ".k")))


map.real = map.real %>% mutate(SNP_knock = ifelse(Knockoff == TRUE, paste0(SNP, ".k"), SNP))
order_in_original <- match(map.real$SNP_knock, map.swapped$SNP_knock)


if(file.exists(paste0(backingfilename_orig, ".bk"))) {
  file.remove(paste0(backingfilename_orig, ".bk"))
}

backingfilename_ordered = paste0(temp_path, "sim_swap_test_ordered", 
                                 "_", my_population, "_res", resolution, "_chr", chromosome)

if(file.exists(paste0(backingfilename_ordered, ".bk"))) {
  file.remove(paste0(backingfilename_ordered, ".bk"))
}

if(file.exists(paste0(backingfilename_ordered, ".rds"))) {
  file.remove(paste0(backingfilename_ordered, ".rds"))
}

G.combined.ordered = big_copy(G.combined, ind.col = order_in_original, backingfile = backingfilename_ordered)
G.combined.ordered <- add_code256(G.combined.ordered, code = obj.bigSNP.real$genotypes$code256)


map.swapped.ordered = map.real.orig %>% 
  mutate(marker.ID = ifelse(endsWith(marker.ID, ".k"), 
                            paste0(gsub(".k", "", marker.ID), ".b"), 
                            paste0(marker.ID, ".a")))

# save new ordered file
# remove files if they existed
if(file.exists((paste0(triple_filename, ".bed")))) {
  file.remove((paste0(triple_filename, ".bed")))
  file.remove((paste0(triple_filename, ".bim")))
  file.remove((paste0(triple_filename, ".fam")))
} 


print("Saving final output structure")
output <- structure(list(genotypes = G.combined.ordered, 
                         fam = obj.bigSNP.real$fam,
                         map = map.swapped.ordered),
                    class = "bigSNP")

print(triple_filename)

snp_writeBed(output, paste0(triple_filename, ".bed"))

# save swapping matrix 
saveRDS(V.swap.matrix, paste0(triple_filename, "_vswap.rds"))

