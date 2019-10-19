import numpy as np
from collections import OrderedDict
from operator import itemgetter

from active_work.read import Dat
from active_work.scde import PDF
from active_work.maths import Histogram, mean_sterr, linspace, logspace,\
    meanStdCut

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

    def nWork(self, n, int_max=None, sum=False):
        """
        Returns active work sums averaged (or not) on packs of `n'.
        """

        workAvegared = []
        for i in self._time0(n, int_max=int_max):
            if sum: workAvegared += [self.activeWork[i:i + n].sum()]
            else: workAvegared += [self.activeWork[i:i + n].mean()]

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

    def nWorkGauss(self, n, *x, cut=None,
        rescaled_to_max=False):
        """
        Returns values of the Gaussian function corresponding to the mean and
        variance of self.nWork(n).
        """

        workSums = self.nWork(n)

        mean, std = meanStdCut(workSums, cut)

        if rescaled_to_max: norm = 1
        else: norm = np.sqrt(2*np.pi*(std**2))

        gauss = lambda y: (
            np.exp(-((y - mean)**2)/(2*(std**2)))
            /norm)

        return np.array(list(map(gauss, x)))

    def corWorkWork(self, tau0=1, n_max=100, int_max=None, max=None, sum=False,
        log=True):
        """

        """

        cor = []
        for n in self._n(n_max=n_max, min=tau0, max=max, log=log):
            if sum == None: worksTot = (lambda l: l - np.mean(l))(
                itemgetter(*(self._time0(n, int_max=int_max) + n))(
                    self.activeWork))
            else: worksTot = (lambda l: l - np.mean(l))(
                self.nWork(n, int_max=int_max, sum=sum))                             # fluctuations of the active wok on intervals of size n
            worksIni = (lambda l: l - np.mean(l))(
                list(map(
                    lambda t: (self.activeWork[t:t + tau0].sum() if sum
                        else self.activeWork[t:t + tau0].mean()),    # fluctuations of the active work averaged on tau0 at the beginning of these intervals
                    self._time0(n, int_max=int_max))))
            # worksIni = (lambda l: l - np.mean(l))(  # fluctuations of the active work at the beginning of these intervals
            #     itemgetter(*self._time0(n, int_max=int_max))(self.activeWork))
            workWork = worksTot*worksIni
            cor += [[n, *mean_sterr(workWork)]]

        return np.array(cor)

    def corWorkOrder(self, n_max=100, int_max=None, max=None, log=False):
        """

        """

        cor = []
        for n in self._n(n_max=n_max, max=max, log=log):
            works = (lambda l: l - np.mean(l))(self.nWork(n, int_max=int_max))  # fluctuations of the active wok on intervals of size n
            orders = (lambda l: l - np.mean(l))(np.array(list(map(              # fluctations of the order parameter norm at the beginning of these intervals
                lambda t: self.getOrderParameter(t, norm=True),
                self.framesWork*self._time0(n, int_max=int_max)))))
            workOrder = works*orders
            cor += [[n, *mean_sterr(workOrder)]]

        return np.array(cor)

    def corOrderOrder(self, n_max=100, int_max=100, max=None, norm=False,
        log=False):
        """

        """

        if max == None: max = self.frames - 1

        if log: space = logspace
        else: space = linspace

        cor = []
        for tau in space(1, max, n_max):
            time0 = list(OrderedDict.fromkeys(np.linspace(
                self.skip*self.framesWork, self.frames - tau - 1, int_max,
                endpoint=True, dtype=int)))
            ordersIni = (lambda l: np.array(l) - np.mean(l, axis=0))(list(map(
                lambda t: self.getOrderParameter(t, norm=norm), time0)))
            ordersFin = (lambda l: np.array(l) - np.mean(l, axis=0))(list(map(
                lambda t: self.getOrderParameter(t + tau, norm=norm), time0)))
            orderOrder = list(map(lambda x, y: np.dot(x, y),
                *(ordersIni, ordersFin)))
            # return ordersIni, ordersFin, orderOrder
            cor += [[tau, *mean_sterr(orderOrder)]]

        return np.array(cor)

    def getWorks(self, tau, n_max=100, init=None):
        """
        Returns array of normalised active works for periods of `tau' frames,
        with a maximum of `n_max', dumping `init' initial frames.
        """

        if init == None: init = int(self.frames/2)

        time0 = np.array(list(OrderedDict.fromkeys(
            np.linspace(init, self.frames - tau - 1, n_max,
                endpoint=True, dtype=int))))

        return np.array(list(map(
            lambda t0: self.getWork(t0, t0 + int(tau)),
            time0)))

    def _time0(self, n, int_max=None):
        """
        Returns list of initial times to coarse-grain the list of active work
        sums in `n' packs.
        """

        time0 = np.linspace(
            self.skip, self.numberWork, int((self.numberWork - self.skip)//n),
            endpoint=False, dtype=int)
        if int_max == None: return time0
        indexes = list(OrderedDict.fromkeys(
            np.linspace(0, time0.size, int_max, endpoint=False, dtype=int)))
        return np.array(itemgetter(*indexes)(time0), ndmin=1)

    def _n(self, n_max=100, min=None, max=None, log=False):
        """
        Returns integers linearly or logarithmically scaled between `min' or 1
        and `max' or int((self.numberWork - self.skip)/2) with `n_max' maximum
        of them.
        """

        if max == None: max = int((self.numberWork - self.skip)/2)

        n_max = int(n_max)
        max = int(max)
        if min == None: min = 1
        min = int(min)

        if log: space = logspace
        else: space = linspace

        return space(min, max, n_max)
