import struct
import numpy as np
import os

from active_work.maths import relative_positions

class Dat:
    """
    Read data files from C library.
    """

    def __init__(self, filename):
        """
        Get data from header.
        """

        self.filename = filename
        self.file = open(self.filename, 'rb')
        self.fileSize = os.path.getsize(filename)

        self.N = self._read('i')            # number of particles
        self.lp = self._read('d')           # persistence length
        self.phi = self._read('d')          # packing fraction
        self.L = self._read('d')            # system size
        self.seed = self._read('i')         # random seed
        self.dt = self._read('d')           # time step
        self.framesWork = self._read('i')   # number of frames on which to sum the active work before dumping

        self.headerLength = self.file.tell()            # length of header in bytes
        self.particleLength = 3*self._bpe('d')          # length the data of a single particle takes in a frame
        self.frameLength = self.N*self.particleLength   # length the data of a single frame takes in a file

        self.numberWork = (self.fileSize - self.headerLength)//(
            self.framesWork*self.frameLength + self._bpe('d'))  # number of cumputed work sums
        self.frames = (self.fileSize - self.headerLength
            - self.numberWork*self._bpe('d'))//self.frameLength # number of frames which the file contains

        if self.fileSize != (
            self.headerLength                   # header
            + self.frames*self.frameLength      # frames
            + self.numberWork*self._bpe('d')):  # work sums
            raise ValueError("Invalid file size.")

    def __del__(self):
        self.file.close()

    def getActiveWork(self):
        """
        Returns array of computed active work sums.
        """

        try: return self.work
        except AttributeError:

            self.work = np.empty(self.numberWork)
            for i in range(self.numberWork):
                self.file.seek(
                    self.headerLength                           # header
                    + self.frameLength                          # frame with index 0
                    + (1 + i)*self.framesWork*self.frameLength  # all following packs of self.framesWork frames
                    + i*self._bpe('d'))                         # previous values of the active work
                self.work[i] = self._read('d')

            return self.work

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
