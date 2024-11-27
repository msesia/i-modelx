#!/bin/bash
#
#SBATCH --job-name=mf4
#SBATCH --partition=mypartition
#SBATCH --mem-per-cpu=200G
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --output=slurm_files/slurm-%j.out

set -e
# Any subsequent(*) commands which fail will cause the shell script to exit immediately

# List of resolutions "1" "2" "3" "4" "5" "6"
RESOLUTION_LIST=("4") #"2" "3" "4"

# set parameters
FDR=0.1
DFMAX=20000 
NCORES=20

echo "DFMAX "$DFMAX" ..."

PARTITION="mypartition"

# set population, phenotype and data paths
POP="british"

PHENO_NAME="whr"
PHENO_FILE="phenotype_file"
REAL_DATA_FILE="real_data_file"


# Utility scripts
RSCRIPT_TO_RUN="Rscript --vanilla utils/male_female_whr.R"

###########################################################################
# KF separately for males and females #
###########################################################################

start=$(date +%s)

  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Running KF separately for males and females"
  echo "----------------------------------------------------------------------------------------------------"

  for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
    echo ""
    echo "Processing at resolution "$RESOLUTION" ..."
    echo ""

    # Compute test statistics
    $RSCRIPT_TO_RUN $POP $RESOLUTION $PHENO_NAME $FDR $DFMAX $NCORES $PHENO_FILE $REAL_DATA_FILE

  done



end=$(date +%s)
echo "Elapsed Time marginal statistics: $(($end-$start)) seconds"
