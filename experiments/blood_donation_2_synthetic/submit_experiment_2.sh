#!/bin/bash

EXP_NUM=1

if [[ $EXP_NUM == 1 ]]; then
  # Parameters
  N_LIST=(100 200 500 1000 2000 5000 10000 20000 50000 80000)
  A_LIST=(1)
  ##CI_LIST=("quantile" "sharp")
  CI_LIST=("quantile")
  NUM_SEEDS=100
fi


#TARGET_LINES=$((1+$BATCH_SIZE*4))

# Slurm parameters
MEMO=2G                             # Memory required (2GB)
TIME=03:00:00                       # Time required (3h)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --account=sesia_658 --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS

OUT_DIR="results2"
mkdir -p $OUT_DIR

for ((J=1; J<=5; J++)); do
  for ((SEED=1; SEED<=$NUM_SEEDS; SEED++)); do
    for N in "${N_LIST[@]}"; do
      for A in "${A_LIST[@]}"; do
        for CI in "${CI_LIST[@]}"; do

          if [[ $CI == "sharp" ]]; then
            MODEL_CLASS_LIST=("glmnet.inter" "glmnet.fast" "glmnet")
            TAIL_LIST=("NA")
            K_LIST=(50)
          else
            MODEL_CLASS_LIST=("glmnet")                
            TAIL_LIST=("right") # "left" "right"
            K_LIST=(90) # 10 90
          fi
          
          for MODEL_CLASS in "${MODEL_CLASS_LIST[@]}"; do
            for TAIL in "${TAIL_LIST[@]}"; do
              for K in "${K_LIST[@]}"; do

                JOBN="n"$N"_a"$A"_"$MODEL_CLASS"_"$CI"_t"$TAIL"_k"$K"_seed"$SEED"_j"$J
                
                OUT_FILE=$OUT_DIR"/"$JOBN".txt"

                RUN=1
                if [[ -f $OUT_FILE ]]; then
                  RUN=0
                fi

                if [[ $RUN == 1 ]]; then
                  # Script to be run
                  SCRIPT="experiment_2.sh $N $A $SEED $MODEL_CLASS $CI $TAIL $K $J"

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
