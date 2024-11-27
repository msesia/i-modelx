#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Merge interactions for different covariates


OUT_BASENAME=$1
RESOLUTION=$2
COV_MIN=$3
COV_MAX=$4
NONZERO_BASEPATH=$5
COV_INT_BASEPATH=$6

# List of chromosomes
COV_LIST=$(seq $COV_MIN $COV_MAX)

# Make list of chromosomes to be merged with the first one
MERGE_LIST=$OUT_BASENAME"_mergelist.txt"
rm -f $MERGE_LIST
touch $MERGE_LIST
for COV in ${COV_LIST[@]}; do
echo $COV
    # Basename for the augmented genotype file for this chromosome
    COV_BASENAME=$COV_INT_BASEPATH"_cov"$COV
    echo $COV_BASENAME >> $MERGE_LIST
done

# Merge the augmented data from all chromosomes
plink \
  --bfile $NONZERO_BASEPATH \
  --merge-list $MERGE_LIST \
  --make-bed \
  --memory 5000 \
  --out $OUT_BASENAME  \
  --keep-allele-order \
  --indiv-sort none

rm -f $OUT_BASENAME".log" 
#$MERGE_LIST


