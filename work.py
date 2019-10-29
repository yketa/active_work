"""
Module work provides classes to compute and analyse active work and active work
autocorrelations and correlations with order parameter.
"""

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

        Parameters
        ----------
        filename : string
            Name of input data file.
        skip : int
            Skip the `skip' first computed values of the active work in the
            following calculations. (default: 1)
            NOTE: This can be changed at any time by setting self.skip.
        """

        super().__init__(filename)  # initialise with super class

        self.skip = skip    # skip the `skip' first measurements of the active work in the analysis

    def nWork(self, n, int_max=None):
        """
        Returns normalised rate of active work averaged on packs of size `n' of
        consecutive individual active works.

        NOTE: Individual active work refers to the normalised rate of active
              work on self.dumpPeriod*self.framesWork consecutive frames and
              stored as element of self.activeWork.

        Parameters
        ----------
        n : int
            Size of packs on which to average active work.
        int_max : int or None
            Maximum number of packs consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of packs.
                  int_max cannot exceed the maximum number of nonoverlapping
                  packs.

        Returns
        -------
        workAvegared : float numpy array
            Array of computed active works.
        """

        workAvegared = []
        for i in self._time0(n, int_max=int_max):
            workAvegared += [self.activeWork[i:i + n].mean()]

        return np.array(workAvegared)

    def nWorkPDF(self, n):
        """
        Returns probability density function of normalised  rate ofactive work
        on packs of `n' of consecutive individual active works.

        NOTE: Individual active work refers to the normalised rate of active
              work on self.dumpPeriod*self.framesWork consecutive frames and
              stored as element of self.activeWork.

        Parameters
        ----------
        n : int
            Size of packs on which to average active work.

        Returns
        -------
        axes : numpy array
            Values at which the probability density function is evaluated.
        pdf : float numpy array
            Values of the probability density function.
        """

        pdf = PDF(self.nWork(n))
        return pdf.axes[0], pdf.pdf

    def nWorkHist(self, n, nBins, vmin=None, vmax=None, log=False,
        rescaled_to_max=False):
        """
        Returns histogram with `nBins' bins of normalised rate of active work on packs of `n' of consecutive individual active works.

        NOTE: Individual active work refers to the normalised rate of active
              work on self.dumpPeriod*self.framesWork consecutive frames and
              stored as element of self.activeWork.

        Parameters
        ----------
        n : int
            Size of packs on which to average active work.
        nBins : int
            Number of bins of the histogram.
        vmin : float
            Minimum value of the bins. (default: self.nWork(n).min())
        vmax : float
            Maximum value of the bins. (default: self.nWork(n).max())
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

        Parameters
        ----------
        n : int
            Size of packs on which to average active work.
        x : float
            Values at which to evaluate the Gaussian function.
        cut : float or None
            Width in units of self.nWork(n).std() to consider when computing
            mean and standard deviation. (see active_work.maths.meanStdCut)
            (default: None)
            NOTE: if cut == None, the width is taken to infinity, i.e. no value
                  is excluded.
        rescaled_to_max : bool
            Rescale function by its computed maximum. (default: False)

        Returns
        -------
        gauss : float numpy array
            Values of the Gaussian function at x.
        """

        workSums = self.nWork(n)

        mean, std = meanStdCut(workSums, cut)

        if rescaled_to_max: norm = 1
        else: norm = np.sqrt(2*np.pi*(std**2))

        gauss = lambda y: (
            np.exp(-((y - mean)**2)/(2*(std**2)))
            /norm)

        return np.array(list(map(gauss, x)))


    def corWorkWorkAve(self,
        tau0=1, n_max=100, int_max=None, min=None, max=None, log=True):
        """
        Compute correlations of work averaged over `tau0' at the beginning of
        an interval and work averaged over the whole interval.
        (see https://yketa.github.io/DAMTP_2019_Wiki/#Active%20Brownian%20particles)

        Parameters
        ----------
        tau0 : int
            Number of consecutive individual active works on which to average
            it. (default: 1)
        n_max : int
            Maximum number of values at which to evaluate the correlation.
            (default: 100)
        int_max : int or None
            Maximum number of different intervals to consider in order to
            compute the mean which appears in the correlation expression.
            (default: None)
            NOTE: if int_max == None, then a maximum number of disjoint
                     intervals will be considered.
        min : int or None
            Minimum value at which to compute the correlation. (default: None)
            NOTE: if min == None then this value is passed as `min' to self.n,
                  otherwise the minimum of `tau0' and `min' is taken.
        max : int or None
            Maximum value at which to compute the correlation. (default: None)
            NOTE: this value is passed as `max' to self.n.
        log : bool
            Logarithmically space values at which the correlations are
            computed. (default: True)

        Returns
        -------
        cor : (5, *) numpy array
            Array of:
                (0) value at which the correlation is computed,
                (1) mean of the computed correlation,
                (2) standard error of the computed correlation,
                (3) standard deviation of the active work computed at the
                    beginning of the interval,
                (4) standard deviation of the active work computed over the
                    interval.
        """

        cor = []
        for n in self._n(n_max=n_max, max=max, log=log,
            min=(None if min == None else np.max([min, tau0]))):
            worksTot = (lambda l: l - np.mean(l))(
                self.nWork(n, int_max=int_max))                     # fluctuations of the active wok on intervals of size n
            worksIni = (lambda l: l - np.mean(l))(
                list(map(
                    lambda t: self.activeWork[t:t + tau0].mean(),   # fluctuations of the active work averaged on tau0 at the beginning of these intervals
                    self._time0(n, int_max=int_max))))
            workWork = worksTot*worksIni
            cor += [[n, *mean_sterr(workWork),
                np.std(worksTot), np.std(worksIni)]]

        return np.array(cor)

    def corWorkWorkIns(self,
        tau0=1, n_max=100, int_max=None, min=None, max=None, log=True):
        """
        Compute correlations of work averaged over `tau0' between different
        times.
        (see https://yketa.github.io/DAMTP_2019_Wiki/#Active%20Brownian%20particles)

        Parameters
        ----------
        tau0 : int
            Number of consecutive individual active works on which to average
            it. (default: 1)
        n_max : int
            Maximum number of values at which to evaluate the correlation.
            (default: 100)
        int_max : int or None
            Maximum number of different intervals to consider in order to
            compute the mean which appears in the correlation expression.
            (default: None)
            NOTE: if int_max == None, then a maximum number of disjoint
                  intervals will be considered.
        min : int or None
            Minimum value at which to compute the correlation. (default: None)
            NOTE: if min == None then min = `tau0', otherwise the minimum of
                  `tau0' and `min' is taken.
        max : int or None
            Maximum value at which to compute the correlation. (default: None)
            NOTE: this value is passed as `max' to self.n.
        log : bool
            Logarithmically space values at which the correlations are
            computed. (default: True)

        Returns
        -------
        cor : (5, *) numpy array
            Array of:
                (0) value at which the correlation is computed,
                (1) mean of the computed correlation,
                (2) standard error of the computed correlation,
                (3) standard deviation of the active work computed at the
                    beginning of the interval,
                (4) standard deviation of the active work computed at the end
                    of the interval.
        """

        cor = []
        for n in self._n(n_max=n_max, log=log,
            min=(2*tau0 if min == None else np.max([tau0 + min, 2*tau0])),
            max=(None if max == None else tau0 + max)):
            worksIni = (lambda l: l - np.mean(l))(
                list(map(
                    lambda t: self.activeWork[t:t + tau0].mean(),   # fluctuations of the active work averaged between t0 and t0 + tau0
                    self._time0(n, int_max=int_max))))
            worksFin = (lambda l: l - np.mean(l))(
                list(map(
                    lambda t: self.activeWork[t + n - tau0:t + n].mean(),       # fluctuations of the active work averaged between t0 and t0 + tau0
                    self._time0(n, int_max=int_max))))
            workWork = worksIni*worksFin
            cor += [[n - tau0, *mean_sterr(workWork),
                np.std(worksIni), np.std(worksFin)]]

        return np.array(cor)

    def corWorkWorkInsBruteForce(self,
        tau0, n_max=100, int_max=None, max=None, log=True):
        """
        Compute correlations of work averaged over `tau0' between different
        times.
        (see https://yketa.github.io/DAMTP_2019_Wiki/#Active%20Brownian%20particles)

        This algorithm computes the correlations more quickly by averaging over
        successive couples of initial and final values of the active work.
        Results of this function should then be taken with care as some other
        unwanted low-time correlations could be picked.

        Parameters
        ----------
        tau0 : int
            Number of consecutive individual active works on which to average
            it. (default: 1)
        n_max : int
            Maximum number of values at which to evaluate the correlation.
            (default: 100)
        int_max : int or None
            Maximum number of different intervals to consider in order to
            compute the mean which appears in the correlation expression.
            (default: None)
            NOTE: if int_max == None, then a maximum number of intervals will
                  intervals will be considered.
        max : int or None
            Maximum value at which to compute the correlation in units of tau0.
            (default: None)
            NOTE: if max == None, the maximum number of values is computed.
        log : bool
            Logarithmically space values at which the correlations are
            computed. (default: True)

        Returns
        -------
        cor : (3, *) numpy array
            Array of:
                (0) value at which the correlation is computed,
                (1) mean of the computed correlation,
                (2) standard error of the computed correlation.
        """

        if log: space = logspace
        else: space = linspace

        if int_max == None: int_max = (self.frames - self.skip)//tau0
        Nsample = int(np.min([(self.frames - self.skip)//tau0, int_max*tau0]))  # size of the sample of consecutive normalised rates of active work to consider
        activeWork = np.array(list(map(                                         # array of consecutive normalised rate of active work averaged of time tau0
            lambda t: self.activeWork[
                self.skip + t*tau0:self.skip + t*tau0 + tau0].mean(),
            range(Nsample))))
        activeWork -= activeWork.mean()                                     # only considering fluctuations to the mean

        lagTimes = space( # array of lag times considered
            1,
            (Nsample - 1) if max == None else int(np.min([max, Nsample - 1])),
            n_max)

        cor = list(map(
            lambda dt: [
                tau0*dt,
                *mean_sterr(
                    (activeWork*np.roll(activeWork, -dt))[:Nsample - dt])],
            lagTimes))

        return np.array(cor)

    def varWorkFromCorWork(self, tau0, n=100, int_max=None, bruteForce=True):
        """
        Compute variance of the active work from its "instantaneous"
        correlations.
        (see https://yketa.github.io/DAMTP_2019_Wiki/#Active%20Brownian%20particles)

        This function is primarily for consistency testing of
        the correlations functions.

        Parameters
        ----------
        tau0 : int
            Number of consecutive individual active works on which to average
            it, and for which correlations will be computed. (default: 1)
        n : int
            Compute variance for tau = i*tau0 with i in {1, ..., n}.
            (default: 100)
        int_max : int or None
            Maximum number of different intervals to consider in order to
            compute the mean which appears in the correlation expression.
            (default: None)
            NOTE: if int_max == None, then a maximum number of intervals will
                  intervals will be considered (joint if bruteForce else
                  disjoint).
        bruteForce : bool
            Use self.corWorkWorkInsBruteForce rather than self.corWorkWorkIns.
            (default: True)

        Returns
        -------
        var : (3, *) numpy array
            Array of:
                (0) value at which the variance is computed,
                (1) mean of the computed variance,
                (2) standard error of the computed variance.
        """

        if bruteForce: corWorkWorkIns = self.corWorkWorkInsBruteForce
        else: corWorkWorkIns = self.corWorkWorkIns

        if bruteForce:
            cor = self.corWorkWorkInsBruteForce(tau0,
                n_max=n, int_max=int_max, max=n - 1, log=False)
        else:
            cor = self.corWorkWorkIns(tau0,
                n_max=n, int_max=int_max, min=tau0, max=(n - 1)*tau0, log=False)

        var0 = mean_sterr((lambda l: (l - l.mean())**2)
            (self.nWork(tau0, int_max=int_max)))

        var = []
        for n0 in range(1, n + 1):
            var += [[n0*tau0, var0[0]/n0, (var0[1]/n0)**2]]
            for i in range(1, n0):
                var[-1][1] += 2*(n0 - i)*cor[i - 1, 1]/(n0**2)
                var[-1][2] += (2*(n0 - i)*cor[i - 1, 2]/(n0**2))**2
            var[-1][2] = np.sqrt(var[-1][2])

        return np.array(var)

    def corWorkOrder(self, n_max=100, int_max=None, min=None, max=None, log=False):
        """

        """

        cor = []
        for n in self._n(n_max=n_max, min=min, max=max, log=log):
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
        sums in packs of `n'.
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
        if min == None: min = 1

        n_max = int(n_max)
        max = int(max)
        min = int(min)

        if log: space = logspace
        else: space = linspace

        return space(min, max, n_max)
