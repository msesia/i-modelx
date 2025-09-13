#!/bin/bash

A_LIST=(0.4 1.0 2.0)
#A_LIST=(0.4 0.8 1.0 1.2 1.4 1.6 1.8 2.0)
SEED_LIST=$(seq 1 100)

# Slurm parameters
MEMO=2G                             # Memory required (2 GB)
TIME=00-02:00:00                    # Time required (10m)
CORE=1                              # Cores required (1)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS

mkdir -p "../data"
mkdir -p "../data/synthetic"
OUT_DIR="../data/synthetic"

for A in "${A_LIST[@]}"; do
  for SEED in $SEED_LIST; do

    JOBN="donations_"$SEED
    OUT_FILE=$OUT_DIR"/dt_donor_mar042018_synthetic_a"$A"_seed"$SEED".csv"
    COMPLETE=0
    if [[ -f $OUT_FILE ]]; then
      COMPLETE=1
    fi

    if [[ $COMPLETE -eq 0 ]]; then
      # Script to be run
      SCRIPT="generate_donations.sh $A $SEED"
      # Define job name for this chromosome
      OUTF=$LOGS"/"$JOBN".out"
      ERRF=$LOGS"/"$JOBN".err"
      # Assemble slurm order for this job
      ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
      # Print order
      echo $ORD
      # Submit order
      #$ORD
      # Run command now
      #./$SCRIPT
    fi
  done
done
