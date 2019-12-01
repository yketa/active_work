"""
Module cloning launches cloning simulations and provides classes to read output
data files from these simulations.
"""

import numpy as np
from numpy import random

from os import path, devnull
from shutil import rmtree as rmr
from shutil import move

from subprocess import Popen

import pickle

from active_work.read import _Read
from active_work.init import get_env, mkdir
from active_work.exponents import float_to_letters
from active_work.maths import mean_sterr

# FUNCTIONS AND CLASSES

class CloningOutput:
    """
    Read and analyse aggregated data from cloning simulations launched with
    active_work.cloning.
    """

    def __init__(self, filename):
        """
        Get data.

        Parameters
        ----------
        filename : string
            Path to data file.
        """

        with open(filename, 'rb') as input:
            (self.tmax,             # dimensionless time simulated
            self.nc,                # number of clones
            self.nRuns,             # number of different runs
            self.initSim,           # number of initial elementary number of iterations to "randomise"
            self.sValues,           # biasing parameters
            self.seed,              # master random seed of master random seeds
            self.seeds,             # master random seeds
            self.N,                 # number of particles in the system
            self.lp,                # dimensionless persistence length
            self.phi,               # packing fraction
            self._tau,              # elementary number of steps
            self.dt,                # time step
            self.SCGF,              # array of different measurements of the SCGF per value of the biasing parameter
            self.activeWork,        # array of different measurements of the active work per value of the biasing parameter
            self.activeWorkForce,   # array of different measurements of the force part of the active work per value of the biasing parameter
            self.activeWorkOri,     # array of different measurements of the orientation part of the active work per value of the biasing parameter
            self.orderParameter,    # array of different measurements of the order parameter per value of the biasing parameter
            self.walltime           # array of different running time per value of the biasing parameter
            ) = pickle.load(input)

        self.tau = self._tau*self.dt    # dimensionless elementary time

    def meanSterr(self):
        """
        Returns array of mean and standard error of measured data.

        Returns
        -------
        SCGF : (self.sValues.size, 3) float Numpy array
            Scaled cumulant generating function.
        activeWork : (self.sValues.size, 3) float Numpy array
            Normalised rate of active work.
        activeWorkForce : (self.sValues.size, 3) float Numpy array
            Force part of the normalised rate of active work.
        activeWorkOri : (self.sValues.size, 3) float Numpy array
            Orientation part of the normalised rate of active work.
        orderParameter : (self.sValues.size, 3) float Numpy array.
            Order parameter.

        NOTE: (0) Biasing parameter.
              (1) Mean.
              (2) Standard error.
        """

        SCGF = np.empty((self.sValues.size, 3))
        activeWork = np.empty((self.sValues.size, 3))
        activeWorkForce = np.empty((self.sValues.size, 3))
        activeWorkOri = np.empty((self.sValues.size, 3))
        orderParameter = np.empty((self.sValues.size, 3))
        for i in range(self.sValues.size):
            SCGF[i] = [
                self.sValues[i], *mean_sterr(self.SCGF[i])]
            activeWork[i] = [
                self.sValues[i], *mean_sterr(self.activeWork[i])]
            activeWorkForce[i] = [
                self.sValues[i], *mean_sterr(self.activeWorkForce[i])]
            activeWorkOri[i] = [
                self.sValues[i], *mean_sterr(self.activeWorkOri[i])]
            orderParameter[i] = [
                self.sValues[i], *mean_sterr(self.orderParameter[i])]

        return SCGF, activeWork, activeWorkForce, activeWorkOri, orderParameter

class _CloningOutput(_Read):
    """
    Read data from a single cloning simulation.
    """

    def __init__(self, filename):
        """
        Get data.

        Parameters
        ----------
        filename : string
            Path to data file.
        """

        # FILE
        super().__init__(filename)

        # HEADER INFORMATION
        self.tmax = self._read('d')         # dimensionless time simulated
        self.nc = self._read('i')           # number of clones
        self.sValue = self._read('d')       # biasing parameter
        self.seed = self._read('i')         # master random seed
        self.nRuns = self._read('i')        # number of different runs
        self.cloneMethod = self._read('i')  # cloning method
        self.initSim = self._read('i')      # number of initial elementary number of iterations to "randomise" the systems
        self.N = self._read('i')            # number of particles in the system
        self.lp = self._read('d')           # dimensionless persistence length
        self.phi = self._read('d')          # packing fraction
        self.tau = self._read('i')          # elementary number of steps
        self.dt = self._read('d')           # time step

        # FILE PARTS LENGTHS
        self.headerLength = self.file.tell()    # length of header in bytes
        self.runLength = 6*self._bpe('d')       # length the data of a run takes

        # FILE CORRUPTION CHECK
        if self.fileSize != self.headerLength + self.nRuns*self.runLength:
            raise ValueError("Invalid data file size.")

        # MEASUREMENTS
        self.psi = np.empty((self.nRuns,))              # scaled cumulant generating function
        self.activeWork = np.empty((self.nRuns,))       # normalised rate of active work
        self.activeWorkForce = np.empty((self.nRuns,))  # force part of the normalised rate of active work
        self.activeWorkOri = np.empty((self.nRuns,))    # orientation part of the normalised rate of active work
        self.orderParameter = np.empty((self.nRuns,))   # order parameter
        self.walltime = np.empty((self.nRuns,))         # time taken for each run
        for i in range(self.nRuns):
            self.psi[i] = self._read('d')
            self.activeWork[i] = self._read('d')
            self.activeWorkForce[i] = self._read('d')
            self.activeWorkOri[i] = self._read('d')
            self.orderParameter[i] = self._read('d')
            self.walltime[i] = self._read('d')

def filename(N, phi, lp, nc, launch):
    """
    Name of simulation output directory.

    Parameters
    ----------
    N : int
        Number of particles in the system.
    phi : float
        Packing fraction.
    lp : float
        Dimensionless persistence length.
    nc : int
        Number of clones.
    launch : int
        Launch identifier.

    Returns
    -------
    name : str
        File name.
    """

    return 'N%s_D%s_L%s_NC%s_E%s' % (*map(float_to_letters,
        (N, phi, lp, nc, launch)),)

# DEFAULT PARAMETERS

_tmax = 1                   # default dimensionless time to simulate
_nc = 10                    # default number of clones
_seed = random.randint(1e7) # default master random seed
_nRuns = 1                  # default number of different runs
_initSim = 1                # default number of initial elementary number of iterations to "randomise" the systems

_sMin = -0.1    # default minimum value of the biasing parameter
_sMax = 0.1     # default maximum value of the biasing parameter
_sNum = 10      # default number of values of the biasing parameter

_threads = -1           # [openMP] default number of threads

_N = 100    # default number of particles in the system
_lp = 5     # default dimension persistence length
_phi = 0.65 # default packing fraction

_tau = 100  # default elementary number of steps
_dt = 0.001 # default time step

_launch = 0 # default launch identifier

_N_cell = 100                                                           # number of particles above which simulations should be launched with a cell list
_exec_dir = path.join(path.dirname(path.realpath(__file__)), 'build')   # default executable directory
_exec_name = ['cloning', 'cloning_cell_list']                           # default executable name without and with a cell list

_out_dir = _exec_dir    # default simulation output directory

# SCRIPT

if __name__ == '__main__':

    # VARIABLE DEFINITIONS

    # CLONING PARAMETERS
    tmax = get_env('TMAX', default=_tmax, vartype=float)        # dimensionless time to simulate
    nc = get_env('NC', default=_nc, vartype=int)                # number of clones
    nRuns = get_env('NRUNS', default=_nRuns, vartype=int)       # number of different runs
    initSim = get_env('INITSIM', default=_initSim, vartype=int) # number of initial elementary number of iterations to "randomise" the systems

    # BIASING PARAMETERS
    sMin = get_env('SMIN', default=_sMin, vartype=float)    # minimum value of the biasing parameter
    sMax = get_env('SMAX', default=_sMax, vartype=float)    # maximum value of the biasing parameter
    sNum = get_env('SNUM', default=_sNum, vartype=int)      # number of values of the biasing parameter
    sValues = np.linspace(sMin, sMax, sNum, endpoint=True)  # array of values of the biasing parameter

    # RANDOM SEEDS
    seed = get_env('SEED', default=_seed, vartype=int)  # master random seed of master random seeds
    random.seed(seed)                                   # set seed
    seeds = random.randint(1e7, size=(sNum,))           # master random seeds

    # OPENMP PARAMETERS
    threads = get_env('THREADS', default=_threads, vartype=int) # number of threads

    # PHYSICAL PARAMETERS
    N = get_env('N', default=_N, vartype=int)           # number of particles in the system
    lp = get_env('LP', default=_lp, vartype=float)      # dimensionless persistence length
    phi = get_env('PHI', default=_phi, vartype=float)   # packing fraction

    # SIMULATION PARAMETERS
    tau = get_env('TAU', default=_tau, vartype=int) # elementary number of steps
    dt = get_env('DT', default=_dt, vartype=float)  # time step

    # EXECUTABLE PARAMETERS
    exec_dir = get_env('EXEC_DIR', default=_exec_dir, vartype=str)      # executable directory
    exec_name = get_env('EXEC_NAME', default=_exec_name[N >= _N_cell],  # executable name
        vartype=str)

    # OUTPUT FILES PARAMETERS
    launch = get_env('LAUNCH', default=_launch, vartype=int)    # launch identifier
    out_dir = get_env('OUT_DIR', default=_out_dir, vartype=str) # output directory
    sim_name = filename(N, phi, lp, nc, launch)                 # simulation output name
    sim_dir = path.join(out_dir, sim_name)                      # simulation output directory name
    mkdir(sim_dir, replace=True)
    tmp_dir = path.join(sim_dir, 'tmp')                         # temporary files directory
    mkdir(tmp_dir, replace=True)
    tmp_template = '%010d.cloning.out'                          # template of temporary files
    out_file = path.join(sim_dir, sim_name + '.dat')            # simulation output file name

    # LAUNCH

    with open(devnull, 'wb') as DevNull:
        procs = [
            Popen([path.join(exec_dir, exec_name)],
                stdout=DevNull, stderr=DevNull, env={
                    'TMAX': str(tmax), 'NC': str(nc), 'SVALUE': str(sValues[i]),
                        'SEED': str(seeds[i]), 'NRUNS': str(nRuns),
                        'INITSIM': str(initSim),
                    'THREADS': str(threads),
                    'N': str(N), 'LP': str(lp), 'PHI': str(phi),
                    'TAU': str(tau), 'DT': str(dt),
                    'FILE': path.join(tmp_dir, tmp_template % i)})
            for i in range(sNum)]       # launch computations
        for proc in procs: proc.wait()  # wait for them to finish

    # CREATING OUTPUT

    # LOAD TEMPORARY FILES
    tmp_out = []
    for i in range(sNum):
        tmp_out += [_CloningOutput(path.join(tmp_dir, tmp_template % i))]

    # ARRAYS OF DATA
    psi = np.array(
        [tmp_out[i].psi for i in range(sNum)])
    activeWork = np.array(
        [tmp_out[i].activeWork for i in range(sNum)])
    activeWorkForce = np.array(
        [tmp_out[i].activeWorkForce for i in range(sNum)])
    activeWorkOri = np.array(
        [tmp_out[i].activeWorkOri for i in range(sNum)])
    orderParameter = np.array(
        [tmp_out[i].orderParameter for i in range(sNum)])
    walltime = np.array(
        [tmp_out[i].walltime for i in range(sNum)])

    # OUT
    with open(out_file, 'wb') as output:
        pickle.dump([
            tmax, nc, nRuns, initSim, sValues,
            seed, seeds,
            N, lp, phi,
            tau, dt,
            psi, activeWork, activeWorkForce, activeWorkOri, orderParameter,
                walltime],
            output)

    # CLEAN
    if get_env('CLEAN', default=True, vartype=bool):
        move(out_file, path.join(out_dir, sim_name + '.dat'))   # move output file to output directory
        rmr(sim_dir)                                            # delete simulation directory
