#!/bin/bash
#
#SBATCH --job-name=whr4
#SBATCH --partition=mypartition
#SBATCH --mem-per-cpu=200G
#SBATCH --time=24:00:00
#SBATCH --output=swapb-%j.out
#SBATCH --output=slurm_files/swapb-%j.out


set -e

# List of covariates
COV_MIN=1
COV_MAX=2
COV_LIST=($(seq $COV_MIN $COV_MAX))

# number of max interactions to be considered
NUM_INTER=1

# List of resolutions 
RESOLUTION_LIST=("4")

# set parameters 
FDR=0.1

DFMAX=20000 
echo "DFMAX "$DFMAX" ..."

NCORES=20
TIMELASSO=48
MEMLASSO=150
DFMAXLASSO=10000
PARTITION="mypartition"

# set population, phenotype, covariate and data paths
POP="british"

PHENO_NAME="whr"
PHENO_FILE="phenotype_file"
REAL_DATA_FILE="real_data_file"

COVAR_ID="sex"

SUFFIX=""


# Utility scripts
COMPUTE_STATS="Rscript --vanilla utils/lasso_swapped.R"
CREATE_INTERACTIONS="Rscript --vanilla utils/create_interactions.R"
MERGE_INTERACTIONS="utils/merge_interactions.sh"
BED_TO_FBM="Rscript --vanilla utils/make_FBM.R"
COMPUTE_STATS_INTERACTION="Rscript --vanilla utils/lasso_on_interactions.R"
COMPUTE_LASSO_ENVIRONMENTS="Rscript --vanilla utils/lasso_sep_each_environment.R"
COMPUTE_LASSO_NO_INTER="Rscript --vanilla utils/lasso_real_data_after_prescreen.R"

# Which operations should we perform?
FLAG_COMPUTE_STATS=1
FLAG_CREATE_INTERACTIONS=1
FLAG_MERGE_INTERACTIONS=1
FLAG_STATS_INTERACTIONS=1
FLAG_LASSO_ENVIRONMENTS=1
FLAG_LASSO_NO_INTER=1
  
###########################
# pre-screening #
###########################

start=$(date +%s)

if [[ $FLAG_COMPUTE_STATS == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Pre-screening"
  echo "----------------------------------------------------------------------------------------------------"

  for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
    echo ""
    echo "Processing at resolution "$RESOLUTION" ..."
    echo ""

    # Compute test statistics
    $COMPUTE_STATS $POP $RESOLUTION $PHENO_NAME $DFMAX $NCORES $PHENO_FILE  

  done

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping pre-screening"
  echo "----------------------------------------------------------------------------------------------------"
fi


end=$(date +%s)
echo "Elapsed time pre-screening: $(($end-$start)) seconds"


#######################
# create interactions #
#######################

start=$(date +%s)


if [[ $FLAG_CREATE_INTERACTIONS == 1 ]]; then
  echo "----------------------------------------------------------------------------------------------------"
  echo "Create interactions for each covariate"
  echo "----------------------------------------------------------------------------------------------------"
  
  for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
    for COV in "${COV_LIST[@]}"; do
      echo ""
      echo "Creating interactions for " $POP " covariate "$COV "at resolution "$RESOLUTION "for "$PHENO_NAME
      echo ""

      # Swap genotypes and knockoffs
      $CREATE_INTERACTIONS $POP $RESOLUTION $COV $PHENO_NAME

    done
  done  
  
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping creation of interactions chunks"
  echo "----------------------------------------------------------------------------------------------------"
fi


##########################################
# merge interactions for each covariate #
##########################################

if [[ $FLAG_MERGE_INTERACTIONS == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Merging interactions for each covariate"
  echo "----------------------------------------------------------------------------------------------------"

  for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
    echo ""
    echo "Processing at resolution "$RESOLUTION" ..."
    echo ""

    # Basename for output FBM
    OUT_BASENAME="path_for_output/UKB_swp_"$POP"_res"$RESOLUTION"_"$PHENO_NAME"_covar_interaction_"$COVAR_ID
    NONZERO_BASEPATH="path_for_prescreened_data/UKB_swp_"$POP"_res"$RESOLUTION"_"$PHENO_NAME"_cov_"$COVAR_ID"_ordered_prescreened"
    COV_INT_BASEPATH="path_for_covariate_interactions/covar_interaction_"$POP"_res"$RESOLUTION"_cov_"$COVAR_ID
    # Combine genotypes and knockoffs into bed
    $MERGE_INTERACTIONS $OUT_BASENAME $RESOLUTION $COV_MIN $COV_MAX $NONZERO_BASEPATH $COV_INT_BASEPATH
    # Convert augmented BED to FBM
    $BED_TO_FBM $OUT_BASENAME $OUT_BASENAME
  done

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping conversion of augmented genotypes into FBM"
  echo "----------------------------------------------------------------------------------------------------"
fi


end=$(date +%s)
echo "Elapsed Time create + merge interactions: $(($end-$start)) seconds"


###########################
#  lasso with interactions #
###########################

start=$(date +%s)

if [[ $FLAG_STATS_INTERACTIONS == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Running lasso with interaction"
  echo "----------------------------------------------------------------------------------------------------"

  for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
    echo ""
    echo "Processing at resolution "$RESOLUTION" ..."
    echo ""

    # Compute test statistics
    $COMPUTE_STATS_INTERACTION $POP $RESOLUTION $PHENO_NAME $DFMAX $NCORES $NUM_INTER $PHENO_FILE

  done

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping lasso with interactions"
  echo "----------------------------------------------------------------------------------------------------"
fi


end=$(date +%s)
echo "Elapsed Time lasso with interactions: $(($end-$start)) seconds"


###################################################
#  lasso separately in each partition #
##################################################

start=$(date +%s)

if [[ $FLAG_LASSO_ENVIRONMENTS == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Lasso separately at each partition "
  echo "----------------------------------------------------------------------------------------------------"

  for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
    echo ""
    echo "Processing at resolution "$RESOLUTION" ..."
    echo ""

    echo "Number interactions " $NUM_INTER
    # Compute test statistics
    $COMPUTE_LASSO_ENVIRONMENTS $POP $RESOLUTION $PHENO_NAME $TIMELASSO $MEMLASSO $DFMAXLASSO $NCORES $PARTITION $SUFFIX $NUM_INTER $PHENO_FILE $REAL_DATA_FILE

  done

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Lasso separately in each partition"
  echo "----------------------------------------------------------------------------------------------------"
fi


end=$(date +%s)
echo "Elapsed time lasso each partition: $(($end-$start)) seconds"


###########################################################################
# knockoff filter - without interactions #
###########################################################################

start=$(date +%s)

if [[ $FLAG_LASSO_NO_INTER == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Computing lasso statistics on (un)/screened data (no interactions) + knockoff filter"
  echo "----------------------------------------------------------------------------------------------------"

  for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
    echo ""
    echo "Processing at resolution "$RESOLUTION" ..."
    echo ""

    # Compute test statistics
    $COMPUTE_LASSO_NO_INTER $POP $RESOLUTION $PHENO_NAME $FDR $DFMAX $NCORES $PHENO_FILE $REAL_DATA_FILE

  done

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping lasso wo interactions and knockoff filter"
  echo "----------------------------------------------------------------------------------------------------"
fi


end=$(date +%s)
echo "Elapsed Time marginal statistics: $(($end-$start)) seconds"
