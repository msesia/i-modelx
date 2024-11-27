#!/usr/bin/env Rscript

# UK Biobank GWAS
#
# Class: R script
#
# Reorder swapped genotypes to be in the same order as original genotypes

args <- commandArgs(trailingOnly = TRUE)

my_population <- as.character(args[1])
resolution <- as.numeric(args[2])
swapped_genotypes_path <- as.character(args[3])
original_genotypes_path <- as.character(args[4])
ordered_output_path <- as.character(args[5])

print(my_population)
print(resolution)
print(swapped_genotypes_path)
print(original_genotypes_path)
print(ordered_output_path)

library(bigsnpr)
library(tidyverse)

###### load swapped data #######  

obj.bigSNP.swapped <- snp_attach(paste0(swapped_genotypes_path,".rds"))
G.swapped = obj.bigSNP.swapped$genotypes
map.swapped = obj.bigSNP.swapped$map %>% 
  separate_wider_delim(marker.ID, ".", names = c("marker.ID", "knock_ind")) %>% 
  mutate(marker.ID = ifelse(knock_ind == "b", paste0(marker.ID, ".k"), marker.ID))

####### load original data ####### 

obj.bigSNP.original <- snp_attach(paste0(original_genotypes_path,".rds"))
G.original = obj.bigSNP.original$genotypes
map.original = obj.bigSNP.original$map

#######  bring swapped data into order of original data #######  

order_in_original <- match(map.original$marker.ID, map.swapped$marker.ID)

backingfilename_ordered = paste0("/tmp/order_merged_chr", 
                                 "_", my_population, "_res", resolution)

if(file.exists(paste0(backingfilename_ordered, ".bk"))) {
  file.remove(paste0(backingfilename_ordered, ".bk"))
}

if(file.exists(paste0(backingfilename_ordered, ".rds"))) {
  file.remove(paste0(backingfilename_ordered, ".rds"))
}

G.swapped.ordered = big_copy(G.swapped, ind.col = order_in_original, backingfile = backingfilename_ordered)

# ordered genotypes should now be in order as in the original
# the map is identical to the original one, so only need to relabel marker.ID
map.swapped.ordered = map.original %>% 
  mutate(marker.ID = ifelse(endsWith(marker.ID, ".k"), 
                            paste0(gsub(".k", "", marker.ID), ".b"), 
                            paste0(marker.ID, ".a")))

head(map.swapped.ordered)

# save as bigsnpr object 
# remove files if they existed
if(file.exists((paste0(ordered_output_path, ".bed")))) {
  file.remove((paste0(ordered_output_path, ".bed")))
  file.remove((paste0(ordered_output_path, ".bim")))
  file.remove((paste0(ordered_output_path, ".fam")))
} 

output <- structure(list(genotypes = G.swapped.ordered, 
                         fam = obj.bigSNP.original$fam,
                         map = map.swapped.ordered),
                    class = "bigSNP")

paste0("Now saving at ", paste0(ordered_output_path, ".bed"))

snp_writeBed(output, paste0(ordered_output_path, ".bed"))

print("Done")
