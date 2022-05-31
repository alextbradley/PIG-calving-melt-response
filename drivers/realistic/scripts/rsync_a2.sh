#!/bin/bash

JOBNO=141

# move to run directory
cd ../run

# make target directory
HOMEDIR=$HOMEROOT/rPIG_${JOBNO}/run
ssh $HOMEHOST "mkdir -p $HOMEDIR"

# rsync peripheral files
rsync -avzL hFac* $HOMEHOST:$HOMEDIR
rsync -avzL Depth* $HOMEHOST:$HOMEDIR
rsync -avzL DRF* $HOMEHOST:$HOMEDIR
rsync -avzL DXG* $HOMEHOST:$HOMEDIR
rsync -avzL DYG* $HOMEHOST:$HOMEDIR
rsync -avzL RAC* $HOMEHOST:$HOMEDIR
rsync -avzL RC* $HOMEHOST:$HOMEDIR
rsync -avzL XC* $HOMEHOST:$HOMEDIR
rsync -avzL YC* $HOMEHOST:$HOMEDIR
rsync -avzL stdout_* $HOMEHOST:$HOMEDIR

# for each diagnostic file, make netcdf and rsync
VARS="
      state2D
      stateTheta
      stateSalt
      stateRho
      stateUvel
      stateVvel
      stateWvel
     "

for VAR in $VARS
do
  rm -rf $VAR.nc
  echo 'seconds since 1979-01-01 00:00:00' > file_list
  ls $VAR.*.data >> file_list
  ./mit2nc
  rsync -avzL $VAR.nc $HOMEHOST:$HOMEDIR
  rm -rf $VAR.nc
done

