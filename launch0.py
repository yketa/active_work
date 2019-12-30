"""
Module launch0 launches simulations with all different parameters.
"""

from active_work.exponents import float_to_letters
from active_work.init import get_env

from numpy.random import randint

from os import path
from subprocess import run

# FUNCTIONS AND CLASSES

def filename(N, epsilon, v0, D, Dr, phi, launch):
    """
    Name of simulation output files.

    Parameters
    ----------
    N : int
        Number of particles in the system.
    epsilon : float
        Coefficient parameter of potential.
    v0 : float
        Self-propulsion velocity.
    D : float
        Translational diffusivity.
    Dr : float
        Rotational diffusivity.
    phi : float
        Packing fraction.
    launch : int
        Launch identifier.

    Returns
    -------
    name : str
        File name.
    """

    return 'N%s_F%s_V%s_T%s_R%s_D%s_E%s.dat0' % (*map(float_to_letters,
        (N, epsilon, v0, D, Dr, phi, launch)),)

# DEFAULT VARIABLES

_N = 10                 # default number of particles in the system
_Dr = 1./40.            # default rotational diffusivity
_v0 = 1                 # default self-propulsion velocity
_phi = 0.65             # default packing fraction

_I = 0  # polydispersity index

_seed = randint(1e7)    # default random seed
_dt = 1e-3              # default time step
_Niter = 5e4            # default number of iterations

_Emin = 1   # default minimum energy at which to stop the minimisation

_launch = 0 # default launch identifier

_nWork = 0  # default number of frames on which to sum the active work before dumping (0 => nWork = lp/dt)
_dump = 1   # default boolean to indicate to dump positions and orientations to output file
_period = 1 # default period of dumping of positions and orientations in number of frames

_N_cell = 100                                                           # number of particles above which simulations should be launched with a cell list
_exec_dir = path.join(path.dirname(path.realpath(__file__)), 'build')   # default executable directory
_exec_name = ['simulation0', 'simulation0_cell_list']                   # default executable name without and with a cell list

_out_dir = _exec_dir    # default simulation output directory

# SCRIPT

if __name__ == '__main__':

    # VARIABLE DEFINITIONS

    # SYSTEM PARAMETERS
    N = get_env('N', default=_N, vartype=int)                       # number of particles in the system
    Dr = get_env('DR', default=_Dr, vartype=float)                  # rotational diffusivity
    epsilon = get_env('EPSILON', default=Dr/3., vartype=float)      # coefficient parameter of potential
    v0 = get_env('V0', default=_v0, vartype=float)                  # self-propulsion velocity
    D = get_env('D', default=epsilon, vartype=float)                # translational diffusivity
    phi = get_env('PHI', default=_phi, vartype=float)               # packing fraction
    I = get_env('I', default=_I, vartype=float)                     # polydispersity index

    # SIMULATION PARAMETERS
    seed = get_env('SEED', default=_seed, vartype=int)      # random seed
    dt = get_env('DT', default=_dt, vartype=float)          # time step
    Niter = get_env('NITER', default=_Niter, vartype=int)   # number of iterations

    # FIRE ALGORITHM PARAMETERS
    Emin = get_env('EMIN', default=_Emin, vartype=float)            # minimum energy at which to stop the minimisation
    iterMax = get_env('ITERMAX', default=int(100/dt), vartype=int)  # maximum number of iterations of the algorithm
    dtmin = get_env('DTMIN', default=dt*1e-3, vartype=float)        # minimum time step at which to stop the algorithm
    dt0 = get_env('DT0', default=dt*1e-1, vartype=float)            # initial time step of the algorith
    dtmax = get_env('DTMAX', default=dt*1e1, vartype=float)         # maximum time step of the algorithm

    # NAMING PARAMETERS
    launch = get_env('LAUNCH', default=_launch, vartype=float)  # launch identifier

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
    out_file = filename(N, epsilon, v0, D, Dr, phi, launch)     # simulation output file name

    # LAUNCH

    run(['setsid', path.join(exec_dir, exec_name)], env={
        'N': str(N), 'EPSILON': str(epsilon), 'V0': str(v0), 'D': str(D),
            'DR': str(Dr), 'PHI': str(phi), 'I': str(I),
        'SEED': str(seed),
        'FILE': path.join(out_dir, out_file),
        'DT': str(dt), 'NITER': str(Niter),
        'NWORK': str(nWork),
        'EMIN': str(Emin), 'ITERMAX': str(iterMax), 'DTMIN': str(dtmin),
            'DT0': str(dt0), 'DTMAX': str(dtmax),
        'DUMP': str(dump), 'PERIOD': str(period)})
