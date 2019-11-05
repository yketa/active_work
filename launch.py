"""
Module launch launches simulations.
"""

from active_work.exponents import float_to_letters
from active_work.init import get_env

from numpy.random import randint

from os import path, getcwd
from subprocess import run

# DEFAULT VARIABLES

_N = 10     # default number of particles in the system
_lp = 40    # default dimensionless persistence length
_phi = 0.65 # default packing fraction

_seed = randint(1e7)    # default random seed
_dt = 1e-3              # default time step
_Niter = 5e4            # default number of iterations

_launch = 0 # default launch identifier

_nWork = 0  # default number of frames on which to sum the active work before dumping (0 => nWork = lp/dt)
_dump = 1   # default boolean to indicate to dump positions and orientations to output file
_period = 1 # default period of dumping of positions and orientations in number of frames

_N_cell = 100                                                           # number of particles above which simulations should be launched with a cell list
_exec_dir = path.join(path.dirname(path.realpath(__file__)), 'build')   # default executable directory
_exec_name = ['simulation', 'simulation_cell_list']                     # default executable name without and with a cell list

_out_dir = _exec_dir    # default simulation output directory

# SCRIPT

if __name__ == '__main__':

    # VARIABLE DEFINITIONS

    # SYSTEM PARAMETERS
    N = get_env('N', default=_N, vartype=int)           # number of particles in the system
    lp = get_env('LP', default=_lp, vartype=float)      # dimensionless persistence length
    phi = get_env('PHI', default=_phi, vartype=float)   # packing fraction

    # SIMULATION PARAMETERS
    seed = get_env('SEED', default=_seed, vartype=int)      # random seed
    dt = get_env('DT', default=_dt, vartype=int)            # time step
    Niter = get_env('NITER', default=_Niter, vartype=int)   # number of iterations

    # NAMING PARAMETERS
    launch = get_env('LAUNCH', default=_launch, vartype=int)    # launch identifier

    # OUTPUT PARAMETERS
    nWork = get_env('NWORK', default=_nWork, vartype=int)       # number of frames on which to sum the active work before dumping (0 => nWork = lp/dt)
    dump = get_env('DUMP', default=_dump, vartype=int)          # boolean to indicate to dump positions and orientations to output file
    period = get_env('PERIOD', default=_period, vartype=int)    # period of dumping of positions and orientations in number of frames

    # EXECUTABLE PARAMETERS
    exec_dir = get_env('EXEC_DIR', default=_exec_dir, vartype=str)      # executable directory
    exec_name = get_env('EXEC_NAME', default=_exec_name[N >= _N_cell],  # executable name
        vartype=str)

    # OUTPUT FILE PARAMETERS
    out_dir = get_env('OUT_DIR', default=_out_dir, vartype=str) # simulation output directory

    # NAMING

    out_file = 'N%s_D%s_L%s_E%s.dat' % (*map(float_to_letters,  # simulation output file name
        (N, phi, lp, launch)),)

    # LAUNCH

    run(['setsid', path.join(exec_dir, exec_name)], env={
        'N': str(N), 'LP': str(lp), 'PHI': str(phi), 'SEED': str(seed),
            'FILE': path.join(out_dir, out_file),
        'DT': str(dt), 'NITER': str(Niter),
        'NWORK': str(nWork),
        'DUMP': str(dump), 'PERIOD': str(period)})