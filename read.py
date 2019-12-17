"""
Module read provides classes to read from data files, especially the ones made
by the simulation script.
"""

import struct
import numpy as np
import os
import pickle

from active_work.maths import relative_positions, angle

class _Read:
    """
    Generic class to read output files.
    """

    def __init__(self, filename):
        """
        Load file.

        Parameters
        ----------
        filename : string
            Path to data file.
        """

        # FILE
        self.filename = filename
        self.file = open(self.filename, 'rb')
        self.fileSize = os.path.getsize(filename)

    def __del__(self):
        self.file.close()

    def _bpe(self, type):
        """
        Returns number of bytes corresponding to type.

        Parameters
        ----------
        type : string
            Type of value to read.
        """

        return struct.calcsize(type)

    def _read(self, type):
        """
        Read element from file with type.

        Parameters
        ----------
        type : string
            Type of value to read.
        """

        return struct.unpack(type, self.file.read(self._bpe(type)))[0]

class Dat(_Read):
    """
    Read data files from simulations.
    """

    def __init__(self, filename):
        """
        Get data from header.

        Parameters
        ----------
        filename : string
            Path to data file.
        """

        # FILE
        super().__init__(filename)

        # HEADER INFORMATION
        self.N = self._read('i')                # number of particles
        self.lp = self._read('d')               # persistence length
        self.phi = self._read('d')              # packing fraction
        self.L = self._read('d')                # system size
        self.seed = self._read('i')             # random seed
        self.dt = self._read('d')               # time step
        self.framesWork = self._read('i')       # number of frames on which to sum the active work before dumping
        self.dumpParticles = self._read('b')    # dump positions and orientations to output file
        self.dumpPeriod = self._read('i')       # period of dumping of positions and orientations in number of frames

        # FILE PARTS LENGTHS
        self.headerLength = self.file.tell()                        # length of header in bytes
        self.particleLength = 3*self._bpe('d')*self.dumpParticles   # length the data of a single particle takes in a frame
        self.frameLength = self.N*self.particleLength               # length the data of a single frame takes in a file
        self.workLength = 4*self._bpe('d')                          # length the data of a single work and order parameter dump take in a file

        # ESTIMATION OF NUMBER OF COMPUTED WORK AND ORDER SUMS AND FRAMES
        self.numberWork = (self.fileSize
            - self.headerLength                                     # header
            - self.frameLength                                      # first frame
            )//(
            self.framesWork*self.frameLength
                + self.workLength)                                  # number of cumputed work sums
        self.frames = 0 if not(self.dumpParticles) else (
            self.fileSize - self.headerLength
            - self.numberWork*self.workLength)//self.frameLength    # number of frames which the file contains

        # FILE CORRUPTION CHECK
        if self.fileSize != (
            self.headerLength                   # header
            + self.frames*self.frameLength      # frames
            + self.numberWork*self.workLength): # work sums
            raise ValueError("Invalid data file size.")

        # COMPUTED NORMALISED RATE OF ACTIVE WORK
        self._loadWork()

    def getWork(self, time0, time1):
        """
        Returns normalised active work between frames `time0' and `time1'.

        Parameters
        ----------
        time0 : int
            Initial frame.
        time1 : int
            Final frame.

        Returns
        -------
        work : float
            Normalised rate of active work.
        """

        work = np.sum(list(map(
            lambda t: self._work(t),
            range(int(time0), int(time1)))))
        work /= self.N*self.dt*(time1 - time0)

        return work

    def getPositions(self, time, *particle, **kwargs):
        """
        Returns positions of particles at time.

        Parameters
        ----------
        time : int
            Frame.
        particle : int
            Indexes of particles.
            NOTE: if none is given, then all particles are returned.

        Optional keyword parameters
        ---------------------------
        centre : (2,) array like
            Returns position relative to `centre'.

        Returns
        -------
        positions : (*, 2) float Numpy array
            Positions at `time'.
        """

        if particle == (): particle = range(self.N)

        positions = np.array(list(map(
            lambda index: self._position(time, index),
            particle)))

        if 'centre' in kwargs:
            return relative_positions(positions, kwargs['centre'], self.L)
        return positions

    def getDisplacements(self, time0, time1, *particle, jump=1):
        """
        Returns displacements of particles between `time0' and `time1'.

        Parameters
        ----------
        time0 : int
            Initial frame.
        time1 : int
            Final frame.
        particle : int
            Indexes of particles.
            NOTE: if none is given, then all particles are returned.
        jump : int
            Period in number of frames at which to check if particles have
            crossed any boundary. (default: 1)
            NOTE: `jump' must be chosen so that particles do not move a distance
                  greater than half the box size during this time.

        Returns
        -------
        displacements : (*, 2) float Numpy array
            Displacements between `time0' and `time1'.
        """

        if particle == (): particle = range(self.N)
        time0 = int(time0)
        time1 = int(time1)
        jump = int(jump)

        displacements = -self.getPositions(time0, *particle)

        increments = np.zeros((len(particle), 2))
        positions1 = -displacements.copy()
        for t in list(range(time0, time1, jump)) + [time1 - 1]:
            positions0 = positions1.copy()
            positions1 = self.getPositions(t + 1, *particle)
            increments += (
                ((positions0 - self.L/2)*(positions1 - self.L/2) < 0)   # if position switches "sign"
                *(np.abs(positions0 - self.L/2) > self.L/4)             # (and) if particle is not in the centre of the box
                *np.sign(positions0 - self.L/2)                         # "sign" of position
                *self.L)

        displacements += positions1 + increments
        return displacements

    def getDistancePositions(self, time, particle0, particle1):
        """
        Returns distance between particles with indexes `particle0' and
        `particle1' at time `time' and their respective positions.

        Parameters
        ----------
        time : int
            Index of frame.
        particle0 : int
            Index of first particle.
        particle1 : int
            Index of second particle.

        Returns
        -------
        dist : float
            Distance between particles.
        pos0 : (2,) float numpy array
            Position of particle0.
        pos1 : (2,) float numpy array
            Position of particle1.
        """

        pos0, pos1 = self.getPositions(time, particle0, particle1)
        return np.sqrt(
            self._diffPeriodic(pos0[0], pos1[0])**2
            + self._diffPeriodic(pos0[1], pos1[1])**2), pos0, pos1

    def getDistance(self, time, particle0, particle1):
        """
        Returns distance between particles with indexes `particle0' and
        `particle1' at time `time'.

        Parameters
        ----------
        time : int
            Index of frame.
        particle0 : int
            Index of first particle.
        particle1 : int
            Index of second particle.

        Returns
        -------
        dist : float
            Distance between particles.
        """

        return self.getDistancePositions(time, particle0, particle1)[0]

    def getOrientations(self, time, *particle):
        """
        Returns orientations of particles at time.

        Parameters
        ----------
        time : int
            Frame
        particle : int
            Indexes of particles.
            NOTE: if none is given, then all particles are returned.

        Returns
        -------
        orientations : (*,) float Numpy array
            Orientations at `time'.
        """

        if particle == (): particle = range(self.N)

        return np.array(list(map(
            lambda index: self._orientation(time, index),
            particle)))

    def getDirections(self, time, *particle):
        """
        Returns self-propulsion vector of particles at time.

        Parameters
        ----------
        time : int
            Frame
        particle : int
            Indexes of particles.
            NOTE: if none is given, then all particles are returned.

        Returns
        -------
        orientations : (*, 2) float Numpy array
            Unitary self-propulsion vectors at `time'.
        """

        if particle == (): particle = range(self.N)

        return np.array(list(map(
            lambda theta: np.array([np.cos(theta), np.sin(theta)]),
            self.getOrientations(time, *particle))))

    def getOrderParameter(self, time, norm=False):
        """
        Returns order parameter, i.e. mean direction, at time.

        Parameters
        ----------
        time : int
            Frame.
        norm : bool
            Return norm of order parameter. (default: False)

        Returns
        -------
        orderParameter : float if `norm' else (2,) float Numpy array
            Order parameter at `time'.
        """

        orderParameter = np.sum(self.getDirections(time), axis=0)/self.N
        if norm: return np.sqrt(np.sum(orderParameter**2))
        return orderParameter

    def getGlobalPhase(self, time):
        """
        Returns global phase at time `time'.

        Parameters
        ----------
        time : int
            Frame.

        Returns
        -------
        phi : float
            Global phase in radians.
        """

        return angle(*self.getOrderParameter(time, norm=False))

    def getTorqueIntegral0(self, time0, time1):
        """
        Returns normalised zeroth integral in the expression of the modified
        active work for control-feedback modified dynamics from `time0' to
        `time1'.
        (see https://yketa.github.io/DAMTP_2019_Wiki/#ABP%20cloning%20algorithm)

        NOTE: Using Stratonovitch convention.

        Parameters
        ----------
        time0 : int
            Initial frame.
        time1 : int
            Final frame.

        Returns
        -------
        torqueIntegral : float
            Normalised integral.
        """

        time0, time1 = int(time0), int(time1)
        if time0 == time1: return 0

        torqueIntegral = 0
        for time in range(time0, time1):
            torqueIntegral += np.sum(
                (self.getOrderParameter(time, norm=True)
                    *np.sin(
                        self.getOrientations(time)
                            - self.getGlobalPhase(time))
                + self.getOrderParameter(time + 1, norm=True)
                    *np.sin(
                        self.getOrientations(time + 1)
                            - self.getGlobalPhase(time + 1)))
                *(self.getOrientations(time + 1) - self.getOrientations(time))
            )/2
            # torqueIntegral += np.sum(
            #     list(map(
            #         lambda i: np.sum(
            #             np.sin(self.getOrientations(time, i)
            #                 - self.getOrientations(time))
            #             + np.sin(self.getOrientations(time + 1, i)
            #                 - self.getOrientations(time + 1)))
            #             *(self.getOrientations(time + 1, i)
            #                 - self.getOrientations(time, i)),
            #         range(self.N)))
            # )/2

        return torqueIntegral/(self.N*(time1 - time0)*self.dt)

    def getTorqueIntegral1(self, time0, time1):
        """
        Returns normalised first integral in the expression of the modified
        active work for control-feedback modified dynamics from `time0' to
        `time1'.
        (see https://yketa.github.io/DAMTP_2019_Wiki/#ABP%20cloning%20algorithm)

        NOTE: Using Stratonovitch convention.

        Parameters
        ----------
        time0 : int
            Initial frame.
        time1 : int
            Final frame.

        Returns
        -------
        torqueIntegral : float
            Normalised integral.
        """

        time0, time1 = int(time0), int(time1)
        if time0 == time1: return 0

        torqueIntegral = 0
        for time in range(time0, time1):
            torqueIntegral += (
                self.getOrderParameter(time, norm=True)**2
                + self.getOrderParameter(time + 1, norm=True)**2)/2

        return torqueIntegral/(time1-time0)

    def getTorqueIntegral2(self, time0, time1):
        """
        Returns normalised second integral in the expression of the modified
        active work for control-feedback modified dynamics from `time0' to
        `time1'.
        (see https://yketa.github.io/DAMTP_2019_Wiki/#ABP%20cloning%20algorithm)

        NOTE: Using Stratonovitch convention.

        Parameters
        ----------
        time0 : int
            Initial frame.
        time1 : int
            Final frame.

        Returns
        -------
        torqueIntegral : float
            Normalised integral.
        """

        time0, time1 = int(time0), int(time1)
        if time0 == time1: return 0

        torqueIntegral = 0
        for time in range(time0, time1):
            torqueIntegral += (
                (self.getOrderParameter(time, norm=True)**2)
                    *np.sum(np.sin(
                        self.getOrientations(time)
                            - self.getGlobalPhase(time))**2)
                + (self.getOrderParameter(time + 1, norm=True)**2)
                    *np.sum(np.sin(
                        self.getOrientations(time + 1)
                            - self.getGlobalPhase(time + 1))**2))/2

        return torqueIntegral/(self.N*(time1-time0))

    def _loadWork(self):
        """
        Loads work from file self.filename + '.work.pickle' if it exists or
        extract it from self.filename and pickle it to
        self.filename + '.work.pickle.
        """

        # ACTIVE WORK

        try:    # try loading

            with open(self.filename + '.work.pickle', 'rb') as workFile:
                self.activeWork = pickle.load(workFile)
                if self.activeWork.size != self.numberWork:
                    raise ValueError("Invalid active work array size.")

        except (FileNotFoundError, EOFError):   # active work file does not exist or file is empty

            # COMPUTE
            self.activeWork = np.empty(self.numberWork)
            for i in range(self.numberWork):
                self.file.seek(
                    self.headerLength                           # header
                    + self.frameLength                          # frame with index 0
                    + (1 + i)*self.framesWork*self.frameLength  # all following packs of self.framesWork frames
                    + i*self.workLength)                        # previous values of the active work
                self.activeWork[i] = self._read('d')

            # DUMP
            with open(self.filename + '.work.pickle', 'wb') as workFile:
                pickle.dump(self.activeWork, workFile)

        # ACTIVE WORK (FORCE)

        try:    # try loading

            with open(self.filename + '.work.force.pickle', 'rb') as workFile:
                self.activeWorkForce = pickle.load(workFile)
                if self.activeWorkForce.size != self.numberWork:
                    raise ValueError("Invalid active work (force) array size.")

        except (FileNotFoundError, EOFError):   # active work (force) file does not exist or file is empty

            # COMPUTE
            self.activeWorkForce = np.empty(self.numberWork)
            for i in range(self.numberWork):
                self.file.seek(
                    self.headerLength                           # header
                    + self.frameLength                          # frame with index 0
                    + (1 + i)*self.framesWork*self.frameLength  # all following packs of self.framesWork frames
                    + i*self.workLength                         # previous values of the active work
                    + self._bpe('d'))                           # value of active work
                self.activeWorkForce[i] = self._read('d')

            # DUMP
            with open(self.filename + '.work.force.pickle', 'wb') as workFile:
                pickle.dump(self.activeWorkForce, workFile)

        # ACTIVE WORK (ORIENTATION)

        try:    # try loading

            with open(self.filename + '.work.ori.pickle', 'rb') as workFile:
                self.activeWorkOri = pickle.load(workFile)
                if self.activeWorkOri.size != self.numberWork:
                    raise ValueError("Invalid active work (ori) array size.")

        except (FileNotFoundError, EOFError):   # active work (orientation) file does not exist or file is empty

            # COMPUTE
            self.activeWorkOri = np.empty(self.numberWork)
            for i in range(self.numberWork):
                self.file.seek(
                    self.headerLength                           # header
                    + self.frameLength                          # frame with index 0
                    + (1 + i)*self.framesWork*self.frameLength  # all following packs of self.framesWork frames
                    + i*self.workLength                         # previous values of the active work
                    + 2*self._bpe('d'))                         # value of active work and force part of active work
                self.activeWorkOri[i] = self._read('d')

            # DUMP
            with open(self.filename + '.work.ori.pickle', 'wb') as workFile:
                pickle.dump(self.activeWorkOri, workFile)

        # ORDER PARAMETER

        try:    # try loading

            with open(self.filename + '.order.pickle', 'rb') as workFile:
                self.orderParameter = pickle.load(workFile)
                if self.orderParameter.size != self.numberWork:
                    raise ValueError("Invalid order parameter array size.")

        except (FileNotFoundError, EOFError):   # order parameter file does not exist or file is empty

            # COMPUTE
            self.orderParameter = np.empty(self.numberWork)
            for i in range(self.numberWork):
                self.file.seek(
                    self.headerLength                           # header
                    + self.frameLength                          # frame with index 0
                    + (1 + i)*self.framesWork*self.frameLength  # all following packs of self.framesWork frames
                    + i*self.workLength                         # previous values of the active work
                    + 3*self._bpe('d'))                         # values of the different parts of active work
                self.orderParameter[i] = self._read('d')

            # DUMP
            with open(self.filename + '.order.pickle', 'wb') as workFile:
                pickle.dump(self.orderParameter, workFile)

    def _position(self, time, particle):
        """
        Returns array of position of particle at time.

        Parameters
        ----------
        time : int
            Frame.
        particle : int
            Index of particle.

        Returns
        -------
        position : (2,) float Numpy array
            Position of `particle' at `time'.
        """

        self.file.seek(
            self.headerLength                                           # header
            + time*self.frameLength                                     # other frames
            + particle*self.particleLength                              # other particles
            + (np.max([time - 1, 0])//self.framesWork)*self.workLength) # active work sums (taking into account the frame with index 0)
        return np.array([self._read('d'), self._read('d')])

    def _orientation(self, time, particle):
        """
        Returns orientation of particle at time.

        Parameters
        ----------
        time : int
            Frame.
        particle : int
            Index of particle.

        Returns
        -------
        orientation : (2,) float Numpy array
            Orientation of `particle' at `time'.
        """

        self.file.seek(
            self.headerLength                                           # header
            + time*self.frameLength                                     # other frames
            + particle*self.particleLength                              # other particles
            + 2*self._bpe('d')                                          # positions
            + (np.max([time - 1, 0])//self.framesWork)*self.workLength) # active work sums (taking into account the frame with index 0)
        return self._read('d')

    def _work(self, time):
        """
        Returns active work between `time' and `time' + 1.

        Parameters
        ----------
        time : int
            Frame.

        Returns
        -------
        work : float
            Normalised rate of active work between `time' and `time' + 1.
        """

        time = int(time)

        work = np.sum(list(map(         # sum over time
            lambda u, dr: np.dot(u,dr), # sum over particles
            *(self.getDisplacements(time, time + 1),
                self.getDirections(time)
                    + self.getDirections(time + 1)))))/2

        return work

    def _diffPeriodic(self, x0, x1):
        """
        Returns algebraic distance from x0 to x1 taking into account periodic
        boundary conditions.

        Parameters
        ----------
        x0 : float
            Coordinate of first point.
        x1 : float
            Coordinate of second point.

        Returns
        -------
        diff : float
            Algebraic distance from x0 to x1.
        """

        diff = x1 - x0
        if np.abs(diff) <= self.L/2: return diff

        diff = (1 - 2*(diff > 0))*np.min([
            np.abs(x0) + np.abs(self.L - x1), np.abs(self.L - x0) + np.abs(x1)])
        return diff
