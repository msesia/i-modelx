#!/bin/bash

ml r/4.0.0
Rscript --vanilla generate_donations.R $1 $2
