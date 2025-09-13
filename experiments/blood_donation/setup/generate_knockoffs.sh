#!/bin/bash

module load python/3.7.6
ml gcc/8.3.0
python3 generate_knockoffs.py $1
