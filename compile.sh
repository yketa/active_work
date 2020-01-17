#!/bin/bash

# Compile all executables.

# ROTORS
make clean && ROTORS=yes make;

for CL in yes no; do  # with and without cell lists

  # SIMULATIONS
  for S0 in yes no; do  # general ABPs and custom model
    make clean && SIM0=$S0 CELLLIST=$CL make;
  done

  # CLONING
  for CD in {0..3}; do  # different controlled dynamics
    make clean && CLONING=yes CONTROLLED_DYNAMICS=$CD CELLLIST=$CL make;
  done

done

# CLEAN
make clean
