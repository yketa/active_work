import struct
import numpy as np
import os

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

        self.N = struct.unpack('i',
            self.file.read(self._bytes_per_element('i')))[0]    # number of particles
        self.lp = struct.unpack('d',
            self.file.read(self._bytes_per_element('d')))[0]    # persistence length
        self.L = struct.unpack('d',
            self.file.read(self._bytes_per_element('d')))[0]    # system size
        self.seed = struct.unpack('i',
            self.file.read(self._bytes_per_element('i')))[0]    # system size

        self.header_length = 2*(
            self._bytes_per_element('i') + self._bytes_per_element('d'))    # length of header in bytes
        self.frames = (os.path.getsize(filename) - self.header_length)//(
            (1 + self.N*5)*self._bytes_per_element('d'))                    # number of frames which the file contains

    def __del__(self):
        self.file.close()

    def getTimeStep(self, time):
        """
        Returns time step at time.
        """

        self.file.seek(
            self.header_length                                  # header
            + time*(1 + 5*self.N)*self._bytes_per_element('d')) # frames

        return struct.unpack('d',
            self.file.read(self._bytes_per_element('d')))[0]

    def getPositions(self, time, *particle, wrapped=True):
        """
        Returns positions of particles at time.
        """

        if particle == (): particle = range(self.N)
        if wrapped: position = self._wrapped_position
        else: position = self._unwrapped_position

        return np.array(list(map(
            lambda index: position(time, index),
            particle)))

    def getOrientations(self, time, *particle):
        """
        Returns orientations if particles at time.
        """

        if particle == (): particle = range(self.N)

        return np.array(list(map(
            lambda index: self._orientation(time, index),
            particle)))

    def _wrapped_position(self, time, particle):
        """
        Returns array of wrapped position of particle at time.
        """

        self.file.seek(
            self.header_length                                  # header
            + time*(1 + 5*self.N)*self._bytes_per_element('d')  # frames
            + self._bytes_per_element('d')                      # time step
            + particle*5*self._bytes_per_element('d'))          # other particles

        return np.array(
            [struct.unpack('d',
                self.file.read(self._bytes_per_element('d')))[0],
            struct.unpack('d',
                self.file.read(self._bytes_per_element('d')))[0]])

    def _unwrapped_position(self, time, particle):
        """
        Returns array of unwrapped position of particle.
        """

        self.file.seek(
            self.header_length                                  # header
            + time*(1 + 5*self.N)*self._bytes_per_element('d')  # frames
            + self._bytes_per_element('d')                      # time step
            + particle*5*self._bytes_per_element('d')           # other particles
            + 2*self._bytes_per_element('d'))                   # wrapped coordinates

        return np.array(
            [struct.unpack('d',
                self.file.read(self._bytes_per_element('d')))[0],
            struct.unpack('d',
                self.file.read(self._bytes_per_element('d')))[0]])

    def _orientation(self, time, particle):
        """
        Returns orientation of particle at time.
        """

        self.file.seek(
            self.header_length                                  # header
            + time*(1 + 5*self.N)*self._bytes_per_element('d')  # frames
            + self._bytes_per_element('d')                      # time step
            + particle*5*self._bytes_per_element('d')           # other particles
            + 4*self._bytes_per_element('d'))                   # position coordinates

        return struct.unpack('d',
            self.file.read(self._bytes_per_element('d')))[0]

    def _bytes_per_element(self, type):
        """
        Returns number of bytes corresponding to type.
        """

        return struct.calcsize(type)
