#!/bin/bash

# Compile all executables.

USE_HEUN=no # use Heun's integrations scheme

# ROTORS
(make clean && ROTORS=yes HEUN=$USE_HEUN make) || exit 0;

for CL in yes no; do  # with and without cell lists

  # SIMULATIONS
  for S0 in yes no; do  # general ABPs and custom model
    (make clean && SIM0=$S0 CELLLIST=$CL HEUN=$USE_HEUN make) || exit 0;
  done

  # CLONING
  for CD in {0..3}; do  # different controlled dynamics
    (make clean && CLONING=yes CONTROLLED_DYNAMICS=$CD CELLLIST=$CL HEUN=$USE_HEUN make) || exit 0;
  done

done

# CLEAN
make clean
