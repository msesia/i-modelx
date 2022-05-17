#!/bin/bash

N_LIST=(80000)
INT_LIST=(0)
#N_LIST=(100 200 500 1000 2000 5000 10000 20000)
SEED_LIST=$(seq 1 100)

# Slurm parameters
MEMO=3G                             # Memory required (3 GB)
TIME=00-01:00:00                    # Time required (1h)
CORE=1                              # Cores required (1)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS

mkdir -p "results/"
OUT_DIR="results"

for N in "${N_LIST[@]}"; do
  for INT in "${INT_LIST[@]}"; do
    for SEED in $SEED_LIST; do

      JOBN="experiment_a"$A"_n"$N"_"$INT"_"$SEED
      OUT_FILE=$OUT_DIR"/summary_real_n"$N"_nint"$INT"_seed"$SEED"_split.csv"
      COMPLETE=0
      if [[ -f $OUT_FILE ]]; then
        COMPLETE=1
      fi

      if [[ $COMPLETE -eq 0 ]]; then
        # Script to be run
        SCRIPT="analysis.sh $N $INT $SEED"
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
