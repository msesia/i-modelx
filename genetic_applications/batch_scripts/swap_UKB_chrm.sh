#!/bin/bash
#
#SBATCH --job-name=wbcm
#SBATCH --partition=mypartition
# Evaluate knockoff goodness of fit
#SBATCH --mem-per-cpu=100G
#SBATCH --time=48:00:00
#SBATCH --output=slurm_files/slurm-%j.out

set -e
# Any subsequent(*) commands which fail will cause the shell script to exit immediately

# List of chromosomes
CHR_MIN=1
CHR_MAX=22
CHR_LIST=($(seq $CHR_MIN $CHR_MAX))

# List of resolutions 
RESOLUTION_LIST=("4")

# set population and data paths
POP="whitenonbritish" 
ORIGINAL_DATA_PATH="original_data_path" 

# Utility scripts
AUGMENT_GENOTYPES="utils/merge_chromosomes.sh"
BED_TO_FBM="Rscript --vanilla utils/make_FBM.R"
ORDER_GENOTYPES="Rscript --vanilla utils/order_chr_merged_file.R"
PLOT_GOF="Rscript --vanilla utils/knockoffs_gof_swapping.R"

# Which operations should we perform?
FLAG_PLOT_GOF=1
FLAG_MAKE_FBM=1
FLAG_ORDER=1


##########################
# Evaluate Knockoffs GOF #
##########################

if [[ $FLAG_PLOT_GOF == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Checking swapped goodness-of-fit"
  echo "----------------------------------------------------------------------------------------------------"

  for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
    for CHR in "${CHR_LIST[@]}"; do
    echo ""
    echo "Computing swapped goodness-of-fit diagnostics for chromosome "$CHR" at resolution "$RESOLUTION" ..."
    echo ""

    STATS_BASENAME="swapped_genotypes/gof/UKB_oct_"$POP"_chr"$CHR"_res"$RESOLUTION
    GROUPS_FILE="group_information/ukb_gen_chr"$CHR"_ibd1_res"$RESOLUTION"_grp.txt"
    OUT_BASENAME="tmp/knock_UKB_gof_plots_swp/knockoffs_swapped_"$POP"chr"$CHR"_res"$RESOLUTION

   echo "----------------------------------------------------------------------------------------------------"
   echo "Computing self similarity diagnostics for chromosome "$CHR" at resolution "$RESOLUTION" ..."
   echo "----------------------------------------------------------------------------------------------------"

    # Compute self-similarity diagnostics
    plink --bfile "swapped_genotypes/chr/UKB_swp_"$POP"_chr"$CHR"_res"$RESOLUTION"_joint_whr" \
          --keep-allele-order \
          --freq \
          --r in-phase --ld-window 2 --ld-window-kb 0 \
          --memory 9000 \
          --out $STATS_BASENAME"_self"

      echo "----------------------------------------------------------------------------------------------------"
      echo "Computing exchangeability diagnostics for chromosome "$CHR" at resolution "$RESOLUTION" ..."
      echo "----------------------------------------------------------------------------------------------------"


    # Compute exchangeability diagnostics
   plink --bfile "swapped_genotypes/chr/UKB_swp_"$POP"_chr"$CHR"_res"$RESOLUTION"_joint_whr" \
          --keep-allele-order \
          --r2 in-phase --ld-window 5000 --ld-window-kb 10000 --ld-window-r2 0.01 \
          --memory 9000 \
          --out $STATS_BASENAME


    echo "----------------------------------------------------------------------------------------------------"
    echo "Make GOF plots "$CHR" at resolution "$RESOLUTION" ..."
    echo "----------------------------------------------------------------------------------------------------"

    # Make GOF plots
    $PLOT_GOF $CHR $RESOLUTION $STATS_BASENAME $GROUPS_FILE $OUT_BASENAME
    done
  done

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping swapped knockoffs goodness-of-fit"
  echo "----------------------------------------------------------------------------------------------------"
fi

#####################
# Merge chromsomes #
####################

if [[ $FLAG_MAKE_FBM == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Converting augmented genotypes into FBM"
  echo "----------------------------------------------------------------------------------------------------"

  for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
    echo ""
    echo "Processing at resolution "$RESOLUTION" ..."
    echo ""

    # Basename for output FBM
    OUT_BASENAME="swapped_genotypes/chr_merged/UKB_swp_"$POP"_res"$RESOLUTION
    CHR_BASEPATH="swapped_genotypes/chr/UKB_swp_"$POP
    CHR_BASENAME_MIN="swapped_genotypes/chr/UKB_swp_"$POP"_chr"$CHR_MIN"_res"$RESOLUTION"_joint_whr"
    GRP_BASEPATH="group_information"
    # Combine genotypes and knockoffs into bed
    $AUGMENT_GENOTYPES $OUT_BASENAME $RESOLUTION $CHR_MIN $CHR_MAX $CHR_BASEPATH $CHR_BASENAME_MIN $GRP_BASEPATH
    # Convert augmented BED to FBM
    $BED_TO_FBM $OUT_BASENAME $OUT_BASENAME
  done

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping conversion of augmented genotypes into FBM"
  echo "----------------------------------------------------------------------------------------------------"
fi

#####################
# Order #
####################

if [[ $FLAG_ORDER == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Ordering swapped genotypes"
  echo "----------------------------------------------------------------------------------------------------"

  for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
    echo ""
    echo "Processing at resolution "$RESOLUTION" ..."
    echo ""
    
    SWAPPED_PATH="swapped_genotypes/chr_merged/UKB_swp_"$POP"_res"$RESOLUTION
    ORDERD_OUTPUT_PATH="swapped_genotypes/chr_merged/UKB_swp_"$POP"_res"$RESOLUTION"_ordered_whr"
    
    $ORDER_GENOTYPES $POP $RESOLUTION $SWAPPED_PATH $ORIGINAL_DATA_PATH $ORDERD_OUTPUT_PATH
  done

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping ordering"
  echo "----------------------------------------------------------------------------------------------------"
fi


