![Polydisperse ABPs with WCA potential](https://github.com/yketa/DAMTP_2019_Wiki/raw/master/Images/header.svg?sanitize=true)

# Active work
### Yann-Edwin Keta, DAMTP, University of Cambridge, 2019

## Introduction

This repository contains scripts, for simulation and analysis purposes, developed for a research project, detailed in **[this wiki](https://yketa.github.io/DAMTP_2019_Wiki)**, concerned with active Brownian particles (ABPs) and the large deviations of the work of self-propelling forces (active work). [[Nemoto *et al.*, *Phys. Rev. E* **99**, 022605 (2019)](https://link.aps.org/doi/10.1103/PhysRevE.99.022605)]

Simulation and cloning scripts are written in C++. Wrapper scripts to launch the latter are written in Python, and other Python classes and functions are available to read and analyse the generated data.

While C++ files can be quite cumbersome, Python wrappers are hopefully more readable and commented enough so that their functioning can be easily understood.

## Requirements

All code was developped and tested on 64-bit linux. C++ cloning scripts necessitate `OpenMP`. Python scripts are written for `python3.*`, import the `active_work` package which necessitates the directory containing this repository to be added to the `$PYTHONPATH`, and rely on the following packages:

- `matplotlib`: plotting,
- `numpy`: mathematical functions and array manipulation,
- `scipy`: grid interpolation ([`scde.py`](https://github.com/yketa/active_work/blob/master/scde.py)),
- `fastkde`: kernel density estimation ([`scde.py`](https://github.com/yketa/active_work/blob/master/scde.py)),

which can be installed by running [`pip.sh`](https://github.com/yketa/active_work/blob/master/pip.sh), provided that `pip` is installed.

## Execution

Compilation of all relevant executables, using `g++`, is possible by running [`compile.sh`](https://github.com/yketa/active_work/blob/master/compile.sh).

Given these have been compiled, they can be executed with the Python scripts listed below.

### Simulations of ABPs

ABP model and simulation procedure is detailed in [this tiddler](https://yketa.github.io/DAMTP_2019_Wiki/#Active%20Brownian%20particles).

- Simulations with custom relations between parameters are launched using [`launch.py`](https://github.com/yketa/active_work/blob/master/launch.py).
- Simulations of general ABPs are launched using [`launch0.py`](https://github.com/yketa/active_work/blob/master/launch0.py).

### Simulations of interacting Brownian rotors

Interacting Brownian rotors model is detailed in [this tiddler](https://yketa.github.io/DAMTP_2019_Wiki/#N-interacting%20Brownian%20rotors).

- Simulations are launched using [`launchR.py`](https://github.com/yketa/active_work/blob/master/launchR.py).

### Cloning of ABPs

Principle and computation scheme of the scaled cumulant generating function (SCGF) of the active work and corresponding averages in the biased ensemble are detailed in [this tiddler](https://yketa.github.io/DAMTP_2019_Wiki/#ABP%20cloning%20algorithm).

- Cloning of trajectories of ABPs systems with custom relations between parameters are launched using [`cloning.py`](https://github.com/yketa/active_work/blob/master/cloning.py).

### Cloning of non-interacting Brownian rotors

Principle and computation scheme of the scaled cumulant generating function (SCGF) of the (squared) polarisation and corresponding averages in the biased ensemble are detailed in [this tiddler](https://yketa.github.io/DAMTP_2019_Wiki/#Brownian%20rotors%20cloning%20algorithm).

- Cloning of trajectories of Brownian rotors are launched using [`cloningR.py`](https://github.com/yketa/active_work/blob/master/cloningR.py).
