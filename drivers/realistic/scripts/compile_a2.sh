#!/bin/bash

# Empty the build directory - but first make sure it exists!
if [ -d "../build" ]; then
  cd ../build
  rm -rf *
else
  echo 'There is no build directory'
  exit 1
fi

# Generate a Makefile
$ROOTDIR/tools/genmake2 -ieee -mpi -mods=../code -of=$M_ROOT/build_options/linux_amd64_gfortran_archer2

# Run the Makefile
make depend
make


