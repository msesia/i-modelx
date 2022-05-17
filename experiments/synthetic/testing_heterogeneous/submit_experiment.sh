#!/bin/bash

# Parameters
N_LIST=(200 300 500 700 1000 1500 2000 3000 5000 10000) # 2000 5000)
P_LIST=(100)
SIGNAL_MEAN_LIST=(400)
SIGNAL_STD_LIST=(0) # 100 200 300)
INTERACTION_LIST=(2)
DELTA_LIST=(5)
CLASS_LIST=("glmnet")

BATCH_LIST=$(seq 1 100)
BATCH_SIZE=1

#TARGET_LINES=$((1+$BATCH_SIZE*4))

# Slurm parameters
MEMO=2G                             # Memory required (2GB)
TIME=01:00:00                       # Time required (1h)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS

OUT_DIR="results"
mkdir -p $OUT_DIR
OUT_DIR="results/testing_heterogeneous"
mkdir -p $OUT_DIR

# Loop over configurations and chromosomes
for N in "${N_LIST[@]}"; do
  for P in "${P_LIST[@]}"; do

    for SIGNAL_MEAN in "${SIGNAL_MEAN_LIST[@]}"; do
      for SIGNAL_STD in "${SIGNAL_STD_LIST[@]}"; do
        for INTERACTION in "${INTERACTION_LIST[@]}"; do
          for DELTA in "${DELTA_LIST[@]}"; do
            for CLASS in "${CLASS_LIST[@]}"; do
              for BATCH in $BATCH_LIST; do

                JOBN="n"$N"_p"$P"_a"$SIGNAL_MEAN"-"$SIGNAL_STD"_i"$INTERACTION"_delta"$DELTA"_linear_"$CLASS"_b"$BATCH
                
                OUT_FILE=$OUT_DIR"/"$JOBN".txt"

                RUN=1
                if [[ -f $OUT_FILE ]]; then
                  RUN=0
                fi
                
                if [[ $RUN == 1 ]]; then
                  # Script to be run
                  SCRIPT="experiment.sh $N $P $SIGNAL_MEAN $SIGNAL_STD $INTERACTION $DELTA $CLASS $BATCH $BATCH_SIZE"
                  # Define job name
                  OUTF=$LOGS"/"$JOBN".out"
                  ERRF=$LOGS"/"$JOBN".err"
                  # Assemble slurm order for this job
                  ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
                  # Print order
                  echo $SCRIPT
                  # Submit order
                  #echo $ORD
                  $ORD
                  # Run command now
                  #./$SCRIPT
                fi

              done
            done
          done
        done
      done
    done
  done
done
