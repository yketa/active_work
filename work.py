import numpy as np
from collections import OrderedDict

from active_work.read import Dat
from active_work.scde import PDF
from active_work.maths import Histogram

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

    def nWork(self, n):
        """
        Returns active work sums averaged on packs of `n'.
        """

        workAvegared = []
        for i in np.linspace(
            self.skip, self.numberWork, int((self.numberWork - self.skip)//n),
            endpoint=False, dtype=int):
            workAvegared += [np.mean(self.activeWork[i:i + n])]

        return np.array(workAvegared)

    def nWorkPDF(self, n):
        """
        Returns probability density function of active work sums averaged on
        packs of `n'.
        """

        pdf = PDF(self.nWork(n))
        return pdf.axes[0], pdf.pdf

    def nWorkHist(self, n, nBins, vmin=None, vmax=None, log=False,
        rescaled_to_max=False):
        """
        Returns histogram with `nBins' bins of active work sums averaged on
        packs of `n'.
        """

        workSums = self.nWork(n)

        if vmin == None: vmin = workSums.min()
        if vmax == None: vmax = workSums.max()
        histogram = Histogram(nBins, vmin, vmax)
        histogram.add_values(*workSums)

        bins = histogram.bins
        hist = histogram.get_histogram()
        if rescaled_to_max: hist /= hist.max()
        if not(log): return bins, hist
        else: return bins[hist > 0], np.log(hist[hist > 0])

    def nWorkGauss(self, n, *x,
        rescaled_to_max=False):
        """
        Returns values of the Gaussian function corresponding to the mean and
        variance of self.nWork(n).
        """

        workSums = self.nWork(n)
        if rescaled_to_max: norm = 1
        else: norm = np.sqrt(2*np.pi*workSums.var())
        gauss = lambda y: (
            np.exp(-((y - workSums.mean())**2)/(2*workSums.var()))
            /norm)

        return np.array(list(map(gauss, x)))

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
