#!/bin/bash

ml r/4.0.0
Rscript --vanilla experiment.R $1 $2 $3
