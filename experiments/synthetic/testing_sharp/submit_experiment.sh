#!/bin/bash

# Parameters
N_LIST=(100 200 500 1000) # 2000 5000)
P_LIST=(20)
SIGNAL_MEAN_LIST=(100)
SIGNAL_STD_LIST=(0) # 100 200 300)
INTERACTION_LIST=(0)
MODEL_LIST=("tree")
CLASS_LIST=("glmnet" "rf")

BATCH_LIST=$(seq 1 10)
BATCH_SIZE=1

#TARGET_LINES=$((1+$BATCH_SIZE*4))

# Slurm parameters
MEMO=2G                             # Memory required (2GB)
TIME=00:30:00                       # Time required (2h)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS

OUT_DIR="results"
mkdir -p $OUT_DIR
OUT_DIR="results/estimation_sharp"
mkdir -p $OUT_DIR

# Loop over configurations and chromosomes
for N in "${N_LIST[@]}"; do
  for P in "${P_LIST[@]}"; do

    for SIGNAL_MEAN in "${SIGNAL_MEAN_LIST[@]}"; do
      for SIGNAL_STD in "${SIGNAL_STD_LIST[@]}"; do
        for INTERACTION in "${INTERACTION_LIST[@]}"; do
          for MODEL in "${MODEL_LIST[@]}"; do
            for CLASS in "${CLASS_LIST[@]}"; do
              for BATCH in $BATCH_LIST; do

                JOBN="n"$N"_p"$P"_a"$SIGNAL_MEAN"-"$SIGNAL_STD"_i"$INTERACTION"_"$MODEL"_"$CLASS"_b"$BATCH
                
                OUT_FILE=$OUT_DIR"/"$JOBN".txt"

                RUN=1
                if [[ -f $OUT_FILE ]]; then
                  RUN=0
                fi
                
                if [[ $RUN == 1 ]]; then
                  # Script to be run
                  SCRIPT="experiment.sh $N $P $SIGNAL_MEAN $SIGNAL_STD $INTERACTION $MODEL $CLASS $BATCH $BATCH_SIZE"
                  # Define job name for this chromosome
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
