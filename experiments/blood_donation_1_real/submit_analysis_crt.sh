#!/bin/bash

N_LIST=(80000)
J_LIST=(1 2 3 4 5)

# Slurm parameters
MEMO=3G                             # Memory required (3 GB)
TIME=00-02:00:00                    # Time required (1h)
CORE=1                              # Cores required (1)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS

mkdir -p "results_crt/"
OUT_DIR="results_crt"

for N in "${N_LIST[@]}"; do
  for J in "${J_LIST[@]}"; do

    JOBN="crt_n"$N"_"$J
    OUT_FILE=$OUT_DIR"/pvalues_real_n"$N"_j"$J".csv"
    COMPLETE=0
    if [[ -f $OUT_FILE ]]; then
      COMPLETE=1
    fi

    if [[ $COMPLETE -eq 0 ]]; then
      # Script to be run
      SCRIPT="analysis_crt.sh $N $J"
      # Define job name for this chromosome
      OUTF=$LOGS"/"$JOBN".out"
      ERRF=$LOGS"/"$JOBN".err"
      # Assemble slurm order for this job
      ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
      # Print order
      echo $ORD
      # Submit order
      $ORD
      # Run command now
      #./$SCRIPT
    fi
  done
done
