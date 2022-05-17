#!/bin/bash

# Parameters
N_LIST=(500 1000 2000)
P_LIST=(10)
SIGNAL_MEAN_LIST=(400)
SIGNAL_STD_LIST=(100)
PROP_TREAT_LIST=(-300 -275 -250 -200 -150 -100 -75 -50 -25 0 25 50 75 100 150 200 250 275 300)
TAIL_LIST=("right")
DELTA_LIST=(200)
#K_LIST=(10 20 30 40 50 60 70 80 90)
K_LIST=(90)

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
OUT_DIR="results/crt_2"
mkdir -p $OUT_DIR

# Loop over configurations and chromosomes
for N in "${N_LIST[@]}"; do
  for P in "${P_LIST[@]}"; do

    for SIGNAL_MEAN in "${SIGNAL_MEAN_LIST[@]}"; do
      for SIGNAL_STD in "${SIGNAL_STD_LIST[@]}"; do
        for PROP_TREAT in "${PROP_TREAT_LIST[@]}"; do
          for TAIL in "${TAIL_LIST[@]}"; do
            for DELTA in "${DELTA_LIST[@]}"; do
              for K in "${K_LIST[@]}"; do

                JOBN="n"$N"_p"$P"_a"$SIGNAL_MEAN"-"$SIGNAL_STD"_t"$PROP_TREAT"_"$TAIL"_"$DELTA"_k"$K
                
                OUT_FILE=$OUT_DIR"/"$JOBN".txt"
                echo $OUT_FILE

                RUN=1
                if [[ -f $OUT_FILE ]]; then
                  RUN=0
                fi
                
                if [[ $RUN == 1 ]]; then
                  # Script to be run
                  SCRIPT="experiment_2.sh $N $P $SIGNAL_MEAN $SIGNAL_STD $PROP_TREAT $TAIL $DELTA $K"
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
