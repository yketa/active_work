#!/bin/bash

# We use this bash script to install the Python packages used in the
# `active_work` packages with the `pip` installer.
# ```
# bash pip.sh
# ```

# Please note that this bash script is used instead of a pip requirements file
# to avoid the installation of fastkde to fail if Cython is not installed first.

# PARAMETERS

PYTHON=python
PIP=pip

# INSTALL

PACKAGES=(    # packages to install (order matters)
  # ALL PURPOSES
  matplotlib  # plotting
  numpy       # mathematical functions and array manipulation
  # SELF-CONSISTENT DENSITY ESTIMATION
  scipy       # grid interpolation
  cython      # has to be installed before fastkde
  fastkde     # kernel density estimation
)
for package in ${PACKAGES[@]}; do
  $PYTHON -m $PIP install $package;
done