#!/bin/bash

SEED_LIST=$(seq 1 100)

# Slurm parameters
MEMO=2G                             # Memory required (2 GB)
TIME=00-00:10:00                    # Time required (10m)
CORE=1                              # Cores required (1)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS

mkdir -p "../data"
mkdir -p "../data/knockoffs"
OUT_DIR="../data/knockoffs"

for SEED in $SEED_LIST; do

  JOBN="knockoffs_"$SEED
  OUT_FILE=$OUT_DIR"/dt_donor_mar042018_simple_N5_seed"$SEED".csv"
  COMPLETE=0
  if [[ -f $OUT_FILE ]]; then
    COMPLETE=1
  fi

  if [[ $COMPLETE -eq 0 ]]; then
    # Script to be run
    SCRIPT="generate_knockoffs.sh $SEED"
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
