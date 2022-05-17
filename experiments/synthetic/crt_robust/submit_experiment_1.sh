#!/bin/bash

EXP_NUM=1

if [[ $EXP_NUM == 1 ]]; then
  # Parameters
  N_LIST=(10 20 50 100 200 500 1000) # 10 20 50 100 200 500 1000 2000 5000
  NUM_BATCHES=10
  SIGNAL_MEAN_LIST=(200)
  SIGNAL_STD_LIST=(100)
  INTERACTION_LIST=(0 1 2)
  CI_LIST=("quantile" "sharp")
fi


#TARGET_LINES=$((1+$BATCH_SIZE*4))

# Slurm parameters
MEMO=2G                             # Memory required (2GB)
TIME=03:00:00                       # Time required (3h)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS

OUT_DIR="results"
mkdir -p $OUT_DIR
OUT_DIR="results/crt_robust_1"
mkdir -p $OUT_DIR

# Loop over configurations and chromosomes
for ((BATCH=1; BATCH<=$NUM_BATCHES; BATCH++)); do
  for N in "${N_LIST[@]}"; do
    for SIGNAL_MEAN in "${SIGNAL_MEAN_LIST[@]}"; do
      for SIGNAL_STD in "${SIGNAL_STD_LIST[@]}"; do
        for INTERACTION in "${INTERACTION_LIST[@]}"; do
          for CI in "${CI_LIST[@]}"; do

            if [[ $CI == "sharp" ]]; then
              MODEL_CLASS_LIST=("glmnet.inter" "glmnet.fast" "glmnet")
              TAIL_LIST=("NA")
              K_LIST=(50)
            else
              MODEL_CLASS_LIST=("glmnet")                
              TAIL_LIST=("left" "right")
              K_LIST=(10 90)
            fi
            
            for MODEL_CLASS in "${MODEL_CLASS_LIST[@]}"; do
              for TAIL in "${TAIL_LIST[@]}"; do
                for K in "${K_LIST[@]}"; do

                  JOBN="n"$N"_p10_a"$SIGNAL_MEAN"-"$SIGNAL_STD"_i"$INTERACTION"_"$MODEL_CLASS"_"$CI"_t"$TAIL"_k"$K"_batch"$BATCH

                  OUT_FILE=$OUT_DIR"/"$JOBN".txt"

                  RUN=1
                  if [[ -f $OUT_FILE ]]; then
                    RUN=0
                  fi

                  if [[ $RUN == 1 ]]; then
                    # Script to be run
                    SCRIPT="experiment_1.sh $N $SIGNAL_MEAN $SIGNAL_STD $INTERACTION $MODEL_CLASS $CI $TAIL $K $BATCH"

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
done
