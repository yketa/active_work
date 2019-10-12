import numpy as np
from collections import OrderedDict

from active_work.read import Dat

class ActiveWork(Dat):
    """
    Compute and analyse active work from simulation data.
    """

    def __init__(self, filename, skip=1):
        """
        Loads file.
        """

        super().__init__(filename)  # initialise with super class

        self.skip = skip    # skip the `skip' first measurements of the active work in the analysis

    def varWork(self, n):
        """
        Computes variance of the active work sums averaged on packs of `n'.
        """

        workAvegared = []
        for i in np.linspace(
            self.skip, self.numberWork, int((self.numberWork - self.skip)//n),
            endpoint=False, dtype=int):
            workAvegared += [np.mean(self.activeWork[i:i + n])]

        return np.var(workAvegared)

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
