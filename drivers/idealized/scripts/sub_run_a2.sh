#!/bin/bash

# set job id
JOBNO=141

# clean run directory and link all required files
./prep_run_a2.sh

# submit the job
sbatch -J APIGi_$JOBNO \
       --account $HECACC \
       run_a2.sh

