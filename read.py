import struct
import numpy as np
import os
from tempfile import TemporaryFile

from active_work.maths import relative_positions

class Dat:
    """
    Read data files from C library.
    """

    def __init__(self, filename):
        """
        Get data from header.
        """

        # FILE
        self.filename = filename
        self.file = open(self.filename, 'rb')
        self.fileSize = os.path.getsize(filename)

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

        # ESTIMATION OF NUMBER OF COMPUTED WORK SUMS AND FRAMES
        self.numberWork = (self.fileSize
            - self.headerLength                                 # header
            - self.frameLength                                  # first frame
            )//(
            self.framesWork*self.frameLength + self._bpe('d'))  # number of cumputed work sums
        self.frames = 0 if not(self.dumpParticles) else (
            self.fileSize - self.headerLength
            - self.numberWork*self._bpe('d'))//self.frameLength # number of frames which the file contains

        # FILE CORRUPTION CHECK
        if self.fileSize != (
            self.headerLength                   # header
            + self.frames*self.frameLength      # frames
            + self.numberWork*self._bpe('d')):  # work sums
            raise ValueError("Invalid file size.")

        # COMPUTED NORMALISED RATE OF ACTIVE WORK
        self.activeWork = np.empty(self.numberWork)
        for i in range(self.numberWork):
            self.file.seek(
                self.headerLength                           # header
                + self.frameLength                          # frame with index 0
                + (1 + i)*self.framesWork*self.frameLength  # all following packs of self.framesWork frames
                + i*self._bpe('d'))                         # previous values of the active work
            self.activeWork[i] = self._read('d')

        # ACTIVE WORK
        self.tempFileActiveWork = TemporaryFile(mode='w+b')
        self.computedActiveWork = []

    def __del__(self):
        self.file.close()
        self.tempFileActiveWork.close()

    def getWork(self, time0, time1):
        """
        Returns normalised active work between frames `time0' and `time1'.
        """

        work = np.sum(list(map(
            lambda t: self._work(t),
            # lambda t: np.sum(list(map(      # sum over time
            #     lambda u, dr: np.dot(u,dr), # sum over particles
            #     *(self.getDisplacements(t, t + 1),
            #         self.getDirections(t) + self.getDirections(t + 1))))),
            range(int(time0), int(time1)))))
        # work /= 2*self.N*self.dt*(time1 - time0)
        work /= self.N*self.dt*(time1 - time0)

        return work

    def getPositions(self, time, *particle, **kwargs):
        """
        Returns positions of particles at time.
        """

        if particle == (): particle = range(self.N)

        positions = np.array(list(map(
            lambda index: self._position(time, index),
            particle)))

        if 'centre' in kwargs:
            return relative_positions(positions, kwargs['centre'], self.L)
        return positions

    def getDisplacements(self, time0, time1, *particle):
        """
        Returns displacements of particles between `time0' and `time1'.
        """

        if particle == (): particle = range(self.N)

        displacements = -self.getPositions(time0, *particle)

        increments = np.zeros((len(particle), 2))
        positions1 = -displacements.copy()
        for t in range(int(time0), int(time1)):
            positions0 = positions1.copy()
            positions1 = self.getPositions(t + 1, *particle)
            increments += (
                ((positions0 - self.L/2)*(positions1 - self.L/2) < 0)   # if position switches "sign"
                *(np.abs(positions0 - self.L/2) > self.L/4)             # (and) if particle is not in the centre of the box
                *np.sign(positions0 - self.L/2)                         # "sign" of position
                *self.L)

        displacements += positions1 + increments
        return displacements

    def getOrientations(self, time, *particle):
        """
        Returns orientations of particles at time.
        """

        if particle == (): particle = range(self.N)

        return np.array(list(map(
            lambda index: self._orientation(time, index),
            particle)))

    def getDirections(self, time, *particle):
        """
        Returns self-propulsion vector of particles at time.
        """

        if particle == (): particle = range(self.N)

        return np.array(list(map(
            lambda theta: np.array([np.cos(theta), np.sin(theta)]),
            self.getOrientations(time, *particle))))

    def getOrderParameter(self, time, norm=False):
        """
        Returns order parameter, i.e. mean direction, at time.
        """

        orderParameter = np.sum(self.getDirections(time), axis=0)/self.N
        if norm: return np.sqrt(np.sum(orderParameter**2))
        return orderParameter

    def _position(self, time, particle):
        """
        Returns array of position of particle at time.
        """

        self.file.seek(
            self.headerLength                                           # header
            + time*self.frameLength                                     # other frames
            + particle*self.particleLength                              # other particles
            + (np.max([time - 1, 0])//self.framesWork)*self._bpe('d'))  # active work sums (taking into account the frame with index 0)
        return np.array([self._read('d'), self._read('d')])

    def _orientation(self, time, particle):
        """
        Returns orientation of particle at time.
        """

        self.file.seek(
            self.headerLength                                           # header
            + time*self.frameLength                                     # other frames
            + particle*self.particleLength                              # other particles
            + 2*self._bpe('d')                                          # positions
            + (np.max([time - 1, 0])//self.framesWork)*self._bpe('d'))  # active work sums (taking into account the frame with index 0)
        return self._read('d')

    def _work(self, time):
        """
        Returns active work between `time' and `time' + 1.
        """

        time = int(time)

        if time in self.computedActiveWork: # active work has already been computed

            self.tempFileActiveWork.seek(self._bpe('d')*
                self.computedActiveWork.index(time))
            return struct.unpack('d',
                self.tempFileActiveWork.read(self._bpe('d')))[0]

        else:   # active work has to be computed and dumped

            work = np.sum(list(map(         # sum over time
                lambda u, dr: np.dot(u,dr), # sum over particles
                *(self.getDisplacements(time, time + 1),
                    self.getDirections(time)
                        + self.getDirections(time + 1)))))/2

            self.tempFileActiveWork.seek(0, 2)                      # go to end of file
            self.tempFileActiveWork.write(struct.pack('d', work))   # dump
            self.computedActiveWork += [time]

            return work

    def _bpe(self, type):
        """
        Returns number of bytes corresponding to type.
        """

        return struct.calcsize(type)

    def _read(self, type):
        """
        Read element from file with type.
        """

        return struct.unpack(type, self.file.read(self._bpe(type)))[0]
