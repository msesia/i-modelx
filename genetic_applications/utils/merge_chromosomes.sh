#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Merge knockoffs for different chromosomes


OUT_BASENAME=$1
RESOLUTION=$2
CHR_MIN=$3
CHR_MAX=$4
CHR_BASEPATH=$5
CHR_BASENAME_MIN=$6
GRP_BASEPATH=$7

# List of chromosomes
CHR_LIST=$(seq $CHR_MIN $CHR_MAX)

echo "CHR LIST IS "$CHR_LIST


echo "Make list of chromosomes to be merged with the first one"
# Make list of chromosomes to be merged with the first one
echo $OUT_BASENAME
MERGE_LIST=$OUT_BASENAME"_mergelist.txt"
echo $MERGE_LIST
rm -f $MERGE_LIST
touch $MERGE_LIST
echo "START LOOP"

for CHR in ${CHR_LIST[@]}; do
  if [ $CHR != $CHR_MIN ]; then
    # Basename for the augmented genotype file for this chromosome
    echo $CHR_BASEPATH
    CHR_BASENAME=$CHR_BASEPATH"_chr"$CHR"_res"$RESOLUTION"_joint_whr" 
    echo $CHR_BASENAME >> $MERGE_LIST
  fi
done

echo "This is merge list "$MERGE_LIST

echo "now merge augmented data from all chromosomes"
# Merge the augmented data from all chromosomes
plink \
  --bfile $CHR_BASENAME_MIN \
  --merge-list $MERGE_LIST \
  --make-bed \
  --memory 5000 \
  --out $OUT_BASENAME  \
  --keep-allele-order \
  --indiv-sort none

rm -f $OUT_BASENAME".log" $MERGE_LIST

echo "make list of snp groups"

# Make list of SNP groups
GRP_LIST=$OUT_BASENAME"_grp.txt"
rm -f $GRP_LIST
touch $GRP_LIST
for CHR in ${CHR_LIST[@]}; do
  # Basename for the augmented genotype file for this chromosome
  GRP_FILE=$GRP_BASEPATH"ukb_gen_chr"$CHR"_ibd1_res"$RESOLUTION"_grp.txt"
  tail -n +2 $GRP_FILE >> $GRP_LIST
done
