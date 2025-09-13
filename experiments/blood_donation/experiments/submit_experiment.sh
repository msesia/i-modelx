#!/bin/bash

A_LIST=(0.4)
#N_LIST=(1000 5000 10000)
N_LIST=(100 200 500 1000 2000 5000 10000 20000 50000 80000)
#N_LIST=(100 200 500 1000 2000 5000 10000 20000)
SEED_LIST=$(seq 1 100)

# Slurm parameters
MEMO=5G                             # Memory required (2 GB)
TIME=00-02:00:00                    # Time required (10m)
CORE=1                              # Cores required (1)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS

mkdir -p "results/"
OUT_DIR="results"

for N in "${N_LIST[@]}"; do
  for A in "${A_LIST[@]}"; do
    for SEED in $SEED_LIST; do

      JOBN="experiment_a"$A"_n"$N"_"$SEED
      OUT_FILE=$OUT_DIR"/summary_synthetic_a"$A"_n"$N"_seed"$SEED"_split.csv"
      COMPLETE=0
      if [[ -f $OUT_FILE ]]; then
        COMPLETE=1
      fi

      if [[ $COMPLETE -eq 0 ]]; then
        # Script to be run
        SCRIPT="experiment.sh $A $N $SEED"
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
done
