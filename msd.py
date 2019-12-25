import numpy as np
from collections import OrderedDict

from active_work.read import Dat, Dat0
from active_work.maths import mean_sterr, logspace

class Msd(Dat):
    """
    Compute and analyse mean square displacement from simulation data.
    """

    def __init__(self, filename, skip=1):
        """
        Load file.

        Parameters
        ----------
        filename : string
            Path to data file.
        skip : int
            Skip the `skip' first frames in the following calculations.
            (default: 1)
            NOTE: This can be changed at any time by setting self.skip.
        """

        super().__init__(filename)  # initialise with super class

        self.skip = skip    # skip the `skip' first frames in the analysis

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
                    displacements[i][0] = self.getDisplacements(    # displacements between time0[i] and dt[0]
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
        Returns value of theoretical mean squared displacement at lag time `dt'.

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

        return 4/(3*self.lp)*dt + 2*self.lp*(
            dt + self.lp*(np.exp(-dt/self.lp) - 1))

class Msd0(Dat0, Msd):
    """
    Compute and analyse mean square displacement from simulation data with all
    different parameters.
    """

    def __init__(self, filename, skip=1):
        """
        Load file.

        Parameters
        ----------
        filename : string
            Path to data file.
        skip : int
            Skip the `skip' first frames in the following calculations.
            (default: 1)
            NOTE: This can be changed at any time by setting self.skip.
        """

        super().__init__(filename)  # initialise with super class

        self.skip = skip    # skip the `skip' first frames in the analysis

    def msd_th(self, dt):
        """
        Returns value of theoretical mean squared displacement at lag time `dt'.

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

        return 4*self.D*dt + (2*(self.v0**2)/self.Dr)*(
            dt + (np.exp(-self.Dr*dt) - 1)/self.Dr)
