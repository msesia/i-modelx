#!/bin/bash

EXP_NUM=1

if [[ $EXP_NUM == 1 ]]; then
  # Parameters
  N_LIST=(10 20 50 100 200 500 1000 2000 5000)
  P_LIST=(10)
  SIGNAL_MEAN_LIST=(400)
  SIGNAL_STD_LIST=(100)
  INTERACTION_LIST=(2)
  DELTA_INT_LIST=(5)
  TAIL_LIST=("right")
  DELTA_LIST=(200)
  K_LIST=(90) #10 20 30 40 50 60 70 80 90)
elif [[ $EXP_NUM == 2 ]]; then
  # Parameters
  N_LIST=(1000)
  P_LIST=(10)
  SIGNAL_MEAN_LIST=(400)
  SIGNAL_STD_LIST=(100)
  INTERACTION_LIST=(2)
  DELTA_INT_LIST=(5)
  TAIL_LIST=("right")
  DELTA_LIST=(0 100 200)
  K_LIST=(10 20 30 40 50 60 70 80 90)
elif [[ $EXP_NUM == 3 ]]; then
  # Parameters
  N_LIST=(100 200 500 1000)
  P_LIST=(10)
  SIGNAL_MEAN_LIST=(400)
  SIGNAL_STD_LIST=(100)
  INTERACTION_LIST=(2)
  DELTA_INT_LIST=(5)
  TAIL_LIST=("right")
  DELTA_LIST=(-200 -100 0 100 200 250 300 400 500 600 700)
  K_LIST=(90)
fi


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
OUT_DIR="results/crt_1"
mkdir -p $OUT_DIR

# Loop over configurations and chromosomes
for N in "${N_LIST[@]}"; do
  for P in "${P_LIST[@]}"; do

    for SIGNAL_MEAN in "${SIGNAL_MEAN_LIST[@]}"; do
      for SIGNAL_STD in "${SIGNAL_STD_LIST[@]}"; do
        for INTERACTION in "${INTERACTION_LIST[@]}"; do
          for DELTA_INT in "${DELTA_INT_LIST[@]}"; do
            for TAIL in "${TAIL_LIST[@]}"; do
              for DELTA in "${DELTA_LIST[@]}"; do
                for K in "${K_LIST[@]}"; do

                  JOBN="n"$N"_p"$P"_a"$SIGNAL_MEAN"-"$SIGNAL_STD"_i"$INTERACTION"_delta"$DELTA_INT"_"$TAIL"_"$DELTA"_k"$K
                  
                  OUT_FILE=$OUT_DIR"/"$JOBN".txt"

                  RUN=1
                  if [[ -f $OUT_FILE ]]; then
                    RUN=0
                  fi
                  
                  if [[ $RUN == 1 ]]; then
                    # Script to be run
                    SCRIPT="experiment_1.sh $N $P $SIGNAL_MEAN $SIGNAL_STD $INTERACTION $DELTA_INT $TAIL $DELTA $K"
                      
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
