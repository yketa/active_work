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

        self.N = struct.unpack('i',
            self.file.read(self._bytes_per_element('i')))[0]    # number of particles
        self.lp = struct.unpack('d',
            self.file.read(self._bytes_per_element('d')))[0]    # persistence length
        self.phi = struct.unpack('d',
            self.file.read(self._bytes_per_element('d')))[0]    # packing fraction
        self.L = struct.unpack('d',
            self.file.read(self._bytes_per_element('d')))[0]    # system size
        self.seed = struct.unpack('i',
            self.file.read(self._bytes_per_element('i')))[0]    # random seed

        self.header_length = (              # length of header in bytes
            self._bytes_per_element('i')    # N
            + self._bytes_per_element('d')  # lp
            + self._bytes_per_element('d')  # phi
            + self._bytes_per_element('d')  # L
            + self._bytes_per_element('i')) # seed
        self.frames = (os.path.getsize(filename) - self.header_length)//(
            (1 + self.N*6)*self._bytes_per_element('d'))                    # number of frames which the file contains

    def __del__(self):
        self.file.close()

    def getTimeStep(self, time):
        """
        Returns time step at time.
        """

        self.file.seek(
            self.header_length                                  # header
            + time*(1 + 6*self.N)*self._bytes_per_element('d')) # frames

        return struct.unpack('d',
            self.file.read(self._bytes_per_element('d')))[0]

    def getActiveWork(self, time, *particle):
        """
        Returns active work of particles between frames `time' - 1 and `time'.
        """

        if particle == (): particle = range(self.N)

        return np.array(list(map(
            lambda index: self._active_work(time, index),
            particle)))

    def getPositions(self, time, *particle, wrapped=True, **kwargs):
        """
        Returns positions of particles at time.
        """

        if particle == (): particle = range(self.N)
        if wrapped: position = self._wrapped_position
        else: position = self._unwrapped_position

        positions = np.array(list(map(
            lambda index: position(time, index),
            particle)))

        if 'centre' in kwargs:
            return relative_positions(positions, kwargs['centre'], self.L)
        return positions

    def getDisplacements(self, time0, time1, *particle):
        """
        Returns displacements between time0 and time1.
        """

        return self.getPositions(time1, *particle, wrapped=False)\
            - self.getPositions(time0, *particle, wrapped=False)

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

    def _active_work(self, time, particle):
        """
        Returns active work of particle between frames `time' - 1 and `time'.
        """

        self.file.seek(
            self.header_length                                  # header
            + time*(1 + 6*self.N)*self._bytes_per_element('d')  # frames
            + self._bytes_per_element('d')                      # time step
            + particle*6*self._bytes_per_element('d'))          # active work

        return struct.unpack('d',
            self.file.read(self._bytes_per_element('d')))[0]

    def _wrapped_position(self, time, particle):
        """
        Returns array of wrapped position of particle at time.
        """

        self.file.seek(
            self.header_length                                  # header
            + time*(1 + 6*self.N)*self._bytes_per_element('d')  # frames
            + self._bytes_per_element('d')                      # time step
            + particle*6*self._bytes_per_element('d')           # other particles
            + self._bytes_per_element('d'))                     # active work


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
            + time*(1 + 6*self.N)*self._bytes_per_element('d')  # frames
            + self._bytes_per_element('d')                      # time step
            + particle*6*self._bytes_per_element('d')           # other particles
            + self._bytes_per_element('d')                      # active work
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
            + time*(1 + 6*self.N)*self._bytes_per_element('d')  # frames
            + self._bytes_per_element('d')                      # time step
            + particle*6*self._bytes_per_element('d')           # other particles
            + self._bytes_per_element('d')                      # active work
            + 4*self._bytes_per_element('d'))                   # position coordinates

        return struct.unpack('d',
            self.file.read(self._bytes_per_element('d')))[0]

    def _bytes_per_element(self, type):
        """
        Returns number of bytes corresponding to type.
        """

        return struct.calcsize(type)

# # FORCES OUPUT
# class DatForces:
#     def __init__(self, filename, N):
#         self.filename = filename
#         self.file = open(self.filename, 'rb')
#         self.N = N
#     def __del__(self):
#         self.file.close()
#     def getForce(self, frame, i, j):
#         frame = int(frame)
#         i, j = sorted([i, j])
#         self.file.seek(struct.calcsize('d')*2*(frame*int((self.N*(self.N - 1)/2)) + int((self.N*(self.N - 1) - (self.N - i)*(self.N - i - 1))/2) + (j - i - 1)))
#         return np.array(
#             [struct.unpack('d',
#                 self.file.read(struct.calcsize('d')))[0],
#             struct.unpack('d',
#                 self.file.read(struct.calcsize('d')))[0]])
#     def getForceSum(self, frame, *index):
#         if index == (): index = range(self.N)
#         force = lambda i: np.sum(np.array(list(map(
#             lambda j: np.sign(j - i)*self.getForce(frame, i, j),
#             np.delete(np.array(range(self.N)), i)))),
#             axis=0)
#         return np.array(list(map(force, index)))
