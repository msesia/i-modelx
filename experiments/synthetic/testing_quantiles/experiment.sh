#!/bin/bash

ml gcc/8.3.0
ml r/4.0.0

Rscript --vanilla experiment.R $1 $2 $3 $4 $5 $6
