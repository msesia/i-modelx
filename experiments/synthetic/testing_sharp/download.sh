#DATA=$1

mkdir -p results

rsync -auv sesia@discovery.usc.edu:/home1/sesia/Workspace/randomized_experiments/code/experiments/estimation_1/results/ results/
