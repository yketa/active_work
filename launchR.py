"""
Module launchR launches simulations of interacting Brownian rotors.
"""

from active_work.exponents import float_to_letters
from active_work.init import get_env

from numpy.random import randint

from os import path
from subprocess import run

# FUNCTIONS AND CLASSES

def filename(N, Dr, g, launch):
    """
    Name of simulation output files.

    Parameters
    ----------
    N : int
        Number of rotors in the system.
    Dr : float
        Rotational diffusivity.
    g : float
        Aligning torque parameter.
    launch : int
        Launch identifier.

    Returns
    -------
    name : str
        File name.
    """

    return 'N%s_R%s_G%s_E%s.datR' % (*map(float_to_letters,
        (N, Dr, g, launch)),)

# DEFAULT VARIABLES

_N = 10     # default number of rotors in the system
_Dr = 1./2. # default rotational diffusivity

_seed = randint(1e7)    # default random seed
_dt = 1e-3              # default time step
_Niter = 5e4            # default number of iterations

_launch = 0 # default launch identifier

_period = 1 # default period of dumping of orientations in number of frames

_exec_dir = path.join(path.dirname(path.realpath(__file__)), 'build')   # default executable directory
_exec_name = 'rotors'                                                   # default executable name

_out_dir = _exec_dir    # default simulation output directory

# SCRIPT

if __name__ == '__main__':

    # VARIABLE DEFINITIONS

    # SYSTEM PARAMETERS
    N = get_env('N', default=_N, vartype=int)       # number of particles in the system
    Dr = get_env('DR', default=_Dr, vartype=float)  # rotational diffusivity
    g = get_env('G', default=-Dr/2, vartype=float)  # aligning torque parameter

    # SIMULATION PARAMETERS
    seed = get_env('SEED', default=_seed, vartype=int)      # random seed
    dt = get_env('DT', default=_dt, vartype=float)          # time step
    Niter = get_env('NITER', default=_Niter, vartype=int)   # number of iterations

    # NAMING PARAMETERS
    launch = get_env('LAUNCH', default=_launch, vartype=float)  # launch identifier

    # OUTPUT PARAMETERS
    period = get_env('PERIOD', default=_period, vartype=int)    # period of dumping of orientations in number of frames

    # EXECUTABLE PARAMETERS
    exec_dir = get_env('EXEC_DIR', default=_exec_dir, vartype=str)      # executable directory
    exec_name = get_env('EXEC_NAME', default=_exec_name, vartype=str)   # executable name

    # OUTPUT FILE PARAMETERS
    out_dir = get_env('OUT_DIR', default=_out_dir, vartype=str) # simulation output directory
    out_file = filename(N, Dr, g, launch)                       # simulation output file name

    # LAUNCH

    run(['setsid', path.join(exec_dir, exec_name)], env={
        'N': str(N), 'DR': str(Dr), 'G': str(g),
        'SEED': str(seed),
        'FILE': path.join(out_dir, out_file),
        'DT': str(dt), 'NITER': str(Niter),
        'PERIOD': str(period)})
