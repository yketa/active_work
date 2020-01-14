"""
Module flow provides classes to compute and analyse displacements and velocities
in order to characterise the flow of systems of ABPs.

(see https://yketa.github.io/DAMTP_2019_Wiki/#ABP%20flow%20characteristics)
"""

import numpy as np
from collections import OrderedDict
from operator import itemgetter

from active_work.read import Dat
from active_work.maths import Distribution, wave_vectors_2D, g2Dto1Dgrid,\
    mean_sterr, logspace

class Displacements(Dat):
    """
    Compute and analyse displacements from simulation data.

    (see https://yketa.github.io/DAMTP_2019_Wiki/#Active%20Brownian%20particles)
    """

    def __init__(self, filename, skip=1):
        """
        Loads file.

        Parameters
        ----------
        filename : string
            Name of input data file.
        skip : int
            Skip the `skip' first computed frames in the following calculations.
            (default: 1)
            NOTE: This can be changed at any time by setting self.skip.
        """

        super().__init__(filename)  # initialise with super class

        self.skip = skip    # skip the `skip' first frames in the analysis

    def nDisplacements(self, dt, int_max=None, jump=1, norm=False):
        """
        Returns array of displacements with lag time `dt'.

        Parameters
        ----------
        dt : int
            Displacement lag time.
        int_max : int or None
            Maximum number of intervals to consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of intervals.
        jump : int
            Period in number of frames at which to check if particles have
            crossed any boundary. (default: 1)
            NOTE: `jump' must be chosen so that particles do not move a distance
                  greater than half the box size during this time.
        norm : bool
            Return norm of displacements rather than 2D displacement.
            (default: False)

        Returns
        -------
        displacements : [not(norm)] (*, self.N, 2) float numpy array
                        [norm] (*, self.N) float numpy array
            Array of computed displacements.
        """

        displacements = []
        for time0 in self._time0(dt, int_max=int_max):
            displacements += [
                self.getDisplacements(time0, time0 + dt, jump=jump, norm=norm)]
        displacements = np.array(displacements)

        return displacements

    def displacementsPDF(self, dt, int_max=None, jump=1):
        """
        Returns probability density function of displacement norm over lag time
        `dt'.

        Parameters
        ----------
        dt : int
            Displacement lag time.
        int_max : int or None
            Maximum number of intervals to consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of intervals.
        jump : int
            Period in number of frames at which to check if particles have
            crossed any boundary. (default: 1)
            NOTE: `jump' must be chosen so that particles do not move a distance
                  greater than half the box size during this time.

        Returns
        -------
        axes : numpy array
            Values at which the probability density function is evaluated.
        pdf : float numpy array
            Values of the probability density function.
        """

        return Distribution(self.nDisplacements(
            dt, int_max=int_max, jump=jump, norm=True)).pdf()

    def displacementsHist(self, dt, nBins, int_max=None, jump=1,
        vmin=None, vmax=None, log=False, rescaled_to_max=False):
        """
        Returns histogram with `nBins' bins of displacement norm over lag time
        `dt'.

        Parameters
        ----------
        dt : int
            Displacement lag time.
        nBins : int
            Number of bins of the histogram.
        int_max : int or None
            Maximum number of intervals to consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of intervals.
        jump : int
            Period in number of frames at which to check if particles have
            crossed any boundary. (default: 1)
            NOTE: `jump' must be chosen so that particles do not move a distance
                  greater than half the box size during this time.
        vmin : float
            Minimum value of the bins. (default: minimum computed displacement)
        vmax : float
            Maximum value of the bins. (default: maximum computed displacement)
        log : bool
            Consider the log of the occupancy of the bins. (default: False)
        rescaled_to_max : bool
            Rescale occupancy of the bins by its maximum over bins.
            (default: False)

        Returns
        -------
        bins : float numpy array
            Values of the bins.
        hist : float numpy array
            Occupancy of the bins.
        """

        return Distribution(self.nDisplacements(
            dt, int_max=int_max, jump=jump, norm=True)).hist(
                nBins, vmin=vmin, vmax=vmax, log=log,
                rescaled_to_max=rescaled_to_max)

    def msd(self, n_max=100, int_max=100, min=None, max=None, jump=1):
        """
        Compute mean square displacement.

        Parameters
        ----------
        n_max : int
            Maximum number of times at which to compute the mean square
            displacement. (default: 100)
        int_max : int
            Maximum number of different intervals to consider when averaging
            the mean square displacement. (default: 100)
        min : int or None
            Minimum time at which to compute the displacement. (default: None)
            NOTE: if min == None, then min = 1.
        max : int or None
            Maximum time at which to compute the displacement. (default: None)
            NOTE: if max == None, then max is taken to be the maximum according
                  to the choice of int_max.
        jump : int
            Compute displacements by treating frames in packs of `jump'. This
            can speed up calculation but also give unaccuracies if
            `jump'*self.dt is of the order of the system size. (default: 1)

        Returns
        -------
        msd_sterr : (3, *) float numpy array
            Array of:
                (0) time at which the mean square displacement is computed,
                (1) mean square displacement,
                (2) standard error of the computed mean square displacement.
        """

        min = 1 if min == None else int(min)
        max = ((self.frames - self.skip - 1)//int_max if max == None
            else int(max))

        dt = logspace(min, max, n_max)                                      # array of lag times
        time0 = np.linspace(self.skip, self.frames - dt.max() - 1, int_max, # array of initial times
            endpoint=False, dtype=int)

        # COMPUTE RELEVANT DISPLACEMENTS FROM DATA
        displacements = np.empty((time0.size, dt.size, self.N, 2))
        for j in range(dt.size):
            if j > 0:
                for i in range(time0.size):
                    displacements[i][j] = (         # displacements between time0[i] and time0[i] + dt[j]
                        displacements[i][j - 1]     # displacements between time0[i] and time0[i] + dt[j - 1]
                        + self.getDisplacements(    # displacements between time0[i] + dt[j - 1] and time0[i] + dt[j]
                            time0[i] + dt[j - 1], time0[i] + dt[j],
                            jump=jump))
            else:
                for i in range(time0.size):
                    displacements[i][0] = self.getDisplacements(    # displacements between time0[i] and time0[i] + dt[0]
                        time0[i], time0[i] + dt[0],
                        jump=jump)

        # COMPUTE MEAN SQUARE DISPLACEMENTS
        msd_sterr = []
        for i in range(dt.size):
            disp = displacements[:, i]
            if self.N > 1:
                disp -= disp.mean(axis=1).reshape(time0.size, 1, 2) # substract mean displacement of particles during each considered intervals
            disp.reshape(time0.size*self.N, 2)
            msd_sterr += [[dt[i], *mean_sterr(np.sum(disp**2, axis=-1))]]

        return np.array(msd_sterr)

    def msd_th(self, dt):
        """
        Returns value of theoretical mean squared displacement at lag time `dt'
        for a single ABP.

        (see https://yketa.github.io/DAMTP_2019_Wiki/#One%20ABP)

        Parameters
        ----------
        dt : float
            Lag time at which to evaluate the theoretical mean squared
            displacement.

        Returns
        -------
        msd : float
            Mean squared displacement.
        """

        if self._isDat0:    # general parameters
            return 4*self.D*dt + (2*(self.v0**2)/self.Dr)*(
                dt + (np.exp(-self.Dr*dt) - 1)/self.Dr)
        else:               # custom relations between parameters
            return 4/(3*self.lp)*dt + 2*self.lp*(
                dt + self.lp*(np.exp(-dt/self.lp) - 1))

    def _time0(self, dt, int_max=None):
        """
        Returns array of initial times to evaluate displacements over lag time
        `dt'.

        Parameters
        ----------
        dt : int
            Displacement lag time.
        int_max : int or None
            Maximum number of initial times to return. (default: None)
            NOTE: if int_max == None, a maximum number of them are returned.

        Returns
        -------
        time0 : (*,) float numpy array
            Array of initial times.
        """

        time0 = np.linspace(
            self.skip, self.frames - 1, int((self.frames - 1 - self.skip)//dt),
            endpoint=False, dtype=int)
        if int_max == None: return time0
        indexes = list(OrderedDict.fromkeys(
            np.linspace(0, time0.size, int_max, endpoint=False, dtype=int)))
        return np.array(itemgetter(*indexes)(time0), ndmin=1)

class Velocities(Dat):
    """
    Compute and analyse velocities from simulation data.

    (see https://yketa.github.io/DAMTP_2019_Wiki/#Active%20Brownian%20particles)
    """

    def __init__(self, filename, skip=1):
        """
        Loads file.

        Parameters
        ----------
        filename : string
            Name of input data file.
        skip : int
            Skip the `skip' first computed frames in the following calculations.
            (default: 1)
            NOTE: This can be changed at any time by setting self.skip.
        """

        super().__init__(filename)  # initialise with super class

        self.skip = skip    # skip the `skip' first frames in the analysis

    def nVelocities(self, int_max=None, norm=False):
        """
        Returns array of velocities.

        Parameters
        ----------
        int_max : int or None
            Maximum number of frames to consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of frames.
                  WARNING: This can be very big.
        norm : bool
            Return norm of velocities rather than 2D velocities.
            (default: False)

        Returns
        -------
        velocities : [not(norm)] (*, self.N, 2) float numpy array
                     [norm] (*, self.N) float numpy array
            Array of computed velocities.
        """

        return np.array(list(map(
            lambda time0: self.getVelocities(time0, norm=norm),
            self._time0(int_max=int_max))))

    def velocitiesPDF(self, int_max=None):
        """
        Returns probability density function velocity norm.

        PARAMETERS
        ----------
        int_max : int or None
            Maximum number of frames to consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of frames.
                  WARNING: This can be very big.

        Returns
        -------
        axes : numpy array
            Values at which the probability density function is evaluated.
        pdf : float numpy array
            Values of the probability density function.
        """

        return Distribution(self.nVelocities(int_max=int_max, norm=True)).pdf()

    def velocitiesHist(self, nBins, int_max=None, vmin=None, vmax=None,
        log=False, rescaled_to_max=False):
        """
        Returns histogram with `nBins' bins of displacement norm over lag time
        `dt'.

        Parameters
        ----------
        nBins : int
            Number of bins of the histogram.
        int_max : int or None
            Maximum number of frames to consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of frames.
                  WARNING: This can be very big.
        vmin : float
            Minimum value of the bins. (default: minimum computed velocity)
        vmax : float
            Maximum value of the bins. (default: maximum computed velocity)
        log : bool
            Consider the log of the occupancy of the bins. (default: False)
        rescaled_to_max : bool
            Rescale occupancy of the bins by its maximum over bins.
            (default: False)

        Returns
        -------
        bins : float numpy array
            Values of the bins.
        hist : float numpy array
            Occupancy of the bins.
        """

        return Distribution(self.nVelocities(int_max=int_max, norm=True)).hist(
            nBins, vmin=vmin, vmax=vmax, log=log,
            rescaled_to_max=rescaled_to_max)

    def energySpectrum(self, int_max=None, nBoxes=None):
        """
        Returns kinetic energy spectrum.

        (see https://yketa.github.io/DAMTP_2019_Wiki/#ABP%20flow%20characteristics)

        Parameters
        ----------
        int_max : int or None
            Maximum number of frames to consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of frames.
                  WARNING: This can be very big.
        nBoxes : int
            Number of boxes in each direction to compute the grid of velocities
            which FFT will be computed. (default: None)
            NOTE: If nBoxes == None, then nBoxes = int(sqrt(self.N)).

        Returns
        -------
        E : (*, 2) float Numpy array
            Array of (k, E(k)).
        """

        if nBoxes == None: nBoxes = int(np.sqrt(self.N))
        nBoxes = int(nBoxes)

        wave_vectors = np.sqrt(np.sum(
            wave_vectors_2D(nBoxes, nBoxes, self.L/nBoxes)**2, axis=-1))    # grid of wave vector norms

        FFTsq = []    # list of squared velocity FFT
        for time, velocity in zip(
            self._time0(int_max=int_max),
            self.nVelocities(int_max=int_max, norm=True)):

            velocitiesFFT = np.fft.fft2(    # FFT of displacement grid
                self.toGrid(time, velocity, nBoxes=nBoxes),
                axes=(0, 1))
            FFTsq += [np.real(np.conj(velocitiesFFT)*velocitiesFFT)]

        return g2Dto1Dgrid(
            wave_vectors*np.mean(FFTsq, axis=0), wave_vectors,
            average_grid=False)

    def _time0(self, int_max=None):
        """
        Returns array of frames at which to compute velocities.

        Parameters
        ----------
        int_max : int or None
            Maximum number of frames to consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of frames.
                  WARNING: This can be very big.

        Returns
        -------
        time0 : (*,) float numpy array
            Array of frames.
        """

        if int_max == None: return np.array(range(self.skip, self.frames - 1))
        return np.linspace(
            self.skip, self.frames - 1, int(int_max), endpoint=False, dtype=int)
