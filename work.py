import numpy as np
from collections import OrderedDict

from active_work.read import Dat

class ActiveWork:
    """
    Compute and analyse active work from simulation data.
    """

    def __init__(self, filename):
        """
        Loads file.
        """

        # FILE
        self.filename = filename
        self.dat = Dat(self.filename)   # data object

    def getWork(self, time0, time1):
        """
        Returns normalised active work between frames `time0' and `time1'.
        """

        return np.sum(
            list(map(
                lambda t: np.sum(self.dat.getActiveWork(t + 1)),    # sum over particles
                range(int(time0), int(time1))))                     # sum over time
            )/(self.dat.N*self.dat.getTimeStep(1)*(time1 - time0))  # normalisation

    def getWorks(self, tau, n_max=100, init=None):
        """
        Returns array of normalised active works for periods of `tau' frames,
        with a maximum of `n_max', dumping `init' initial frames.
        """

        if init == None: init = int(self.dat.frames/2)

        time0 = np.array(list(OrderedDict.fromkeys(
            np.linspace(init, self.dat.frames - tau - 1, n_max,
                endpoint=True, dtype=int))))

        return np.array(list(map(
            lambda t0: self.getWork(t0, t0 + int(tau)),
            time0)))
