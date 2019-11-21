"""
Module workorder provides classes to conjointly and analyse active work and
order parameter.
"""

import numpy as np

from active_work.work import ActiveWork
from active_work.maths import Histogram3D
from active_work.scde import PDF

class WorkOrder(ActiveWork):
    """
    Conjointly compute and analyse active work and order parameter.

    (see active_work.work.ActiveWork)
    (see https://yketa.github.io/DAMTP_2019_Wiki/#ABP%20work%20and%20order%20LDP)
    """

    def __init__(self, filename, workPart='all', skip=1):
        """
        Loads file.

        Parameters
        ----------
        filename : string
            Name of input data file.
        workPart : string
            Part of the active work to consider in computations:
                * 'all': active work,
                * 'force': force part of the active work,
                * 'orientation': orientation part of the active work,
                * 'noise': noise part of the active work.
            (default: 'all')
            NOTE: This can be changed at any time by calling self._setWorkPart.
        skip : int
            Skip the `skip' first computed values of the active work in the
            following calculations. (default: 1)
            NOTE: This can be changed at any time by setting self.skip.
        """

        super().__init__(filename, workPart=workPart, skip=skip)    # initialise with super class

    def nWorkOrder(self, n, int_max=None):
        """
        Returns normalised rate of active work and order parameter averaged on
        packs of size `n' of consecutive individual active works and order
        parameters.

        NOTE: Individual active work refers to the normalised rate of active
              work and order parameter on self.dumpPeriod*self.framesWork
              consecutive frames and stored as element of self.workArray and
              self.orderParameter.

        Parameters
        ----------
        n : int
            Size of packs on which to average active work and order parameter.
        int_max : int or None
            Maximum number of packs consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of packs.
                  int_max cannot exceed the maximum number of nonoverlapping
                  packs.

        Returns
        -------
        workAveraged : float numpy array
            Array of computed active works.
        orderAveraged : float numpy array
            Array of computed order parameters.
        """

        workAveraged = []
        orderAveraged = []
        for i in self._time0(n, int_max=int_max):
            workAveraged += [self.workArray[i:i + n].mean()]
            orderAveraged += [self.orderParameter[i:i + n].mean()]

        return np.array(workAveraged), np.array(orderAveraged)

    def SCGF(self, *s, n=1, int_max=None):
        """
        Returns scaled cumulant generating function from active work averaged on
        packs of size `n' of consecutive individual measures at biasing
        parameter `s'.

        (see https://yketa.github.io/DAMTP_2019_Wiki/#ABP%20work%20and%20order%20LDP)

        Parameters
        ----------
        s : float
            Biasing parameter.
        n : int
            Size of packs on which to average active work.
            (default: 1)
        int_max : int or None
            Maximum number of packs consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of packs.
                  int_max cannot exceed the maximum number of nonoverlapping
                  packs.

        Returns
        -------
        tau : float
            Averaging time in absolute dimensionless units.
        psi : float Numpy array
            Scaled cumulant generating function at `s'.
        """

        workArray = super().nWork(n, int_max=int_max)   # only computation of the work is needed
        tau = n*self.dt*self.dumpPeriod*self.framesWork

        return tau, np.array(list(map(
            lambda _s: np.log(np.mean(np.exp(-_s*tau*workArray)))/tau,
            s)))

    def sWork(self, *s, n=1, int_max=None):
        """
        Returns averaged active work in biased ensemble from active work
        averaged on packs of size `n' of consecutive individual measures at
        biasing parameter `s'.

        (see https://yketa.github.io/DAMTP_2019_Wiki/#ABP%20work%20and%20order%20LDP)

        Parameters
        ----------
        s : float
            Biasing parameter.
        n : int
            Size of packs on which to average active work.
            (default: 1)
        int_max : int or None
            Maximum number of packs consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of packs.
                  int_max cannot exceed the maximum number of nonoverlapping
                  packs.

        Returns
        -------
        tau : float
            Averaging time in absolute dimensionless units.
        work : float Numpy array
            Averaged active work at `s'.
        """

        workArray = super().nWork(n, int_max=int_max)   # only computation of the work is needed
        tau = n*self.dt*self.dumpPeriod*self.framesWork

        return tau, np.array(list(map(
            lambda _s: np.mean(workArray*np.exp(-_s*tau*workArray))/(
                np.mean(np.exp(-_s*tau*workArray))),
            s)))

    def sOrder(self, *s, n=1, int_max=None):
        """
        Returns averaged order parameter in biased ensemble from active work and
        order parameter averaged on packs of size `n' of consecutive individual
        measures at biasing parameter `s'.

        (see https://yketa.github.io/DAMTP_2019_Wiki/#ABP%20work%20and%20order%20LDP)

        Parameters
        ----------
        s : float
            Biasing parameter.
        n : int
            Size of packs on which to average active work and order parameter.
            (default: 1)
        int_max : int or None
            Maximum number of packs consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of packs.
                  int_max cannot exceed the maximum number of nonoverlapping
                  packs.

        Returns
        -------
        tau : float
            Averaging time in absolute dimensionless units.
        order : float Numpy array
            Averaged order parameter at `s'.
        """

        workArray, orderArray = self.nWorkOrder(n, int_max=int_max)
        tau = n*self.dt*self.dumpPeriod*self.framesWork

        return tau, np.array(list(map(
            lambda _s: np.mean(orderArray*np.exp(-_s*tau*workArray))/(
                np.mean(np.exp(-_s*tau*workArray))),
            s)))

    def getHistogram3D(self, Nbins, n=1, int_max=None,
        work_min=None, work_max=None, order_min=None, order_max=None):
        """
        Returns 3D histogram of work and order.

        Parameters
        ----------
        Nbins : int or 2-uple-like of int
            Number of histogram bins for active work and order parameter.
        n : int
            Size of packs on which to average active work and order parameter.
        int_max : int or None
            Maximum number of packs consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of packs.
                  int_max cannot exceed the maximum number of nonoverlapping
                  packs.
        work_min : float or None
            Minimum value for the active work. (default: None)
            NOTE: if work_min == None, then the minimum value in the array is
                  taken.
        work_max : float or None
            Maximum value for the active work. (default: None)
            NOTE: if work_max == None, then the maximum value in the array is
                  taken.
        order_min : float or None
            Minimum value for the order parameter. (default: None)
            NOTE: if order_min == None, then the minimum value in the array is
                  taken.
        order_max : float or None
            Maximum value for the order parameter. (default: None)
            NOTE: if order_max == None, then the maximum value in the array is
                  taken.

        Returns
        -------
        hist : (Nbins.prod(), 3) float Numpy array
            Values of the histogram:
                (0) Active work bin.
                (1) Order parameter bin.
                (2) Proportion.
        """

        workArray, orderParameter = self.nWorkOrder(n, int_max=int_max)

        if work_min == None: work_min = np.min(workArray)
        if work_max == None: work_max = np.max(workArray)
        if order_min == None: order_min = np.min(orderParameter)
        if order_max == None: order_max = np.max(orderParameter)

        histogram = Histogram3D(Nbins,
            (work_min, order_min), (work_max, order_max),
            log=False)
        histogram.values = list(zip(workArray, orderParameter))

        return histogram.get_histogram()

    def getHistogram3DSC(self, n=1, int_max=None):
        """
        Returns 3D histogram computed via self-consistent density estimation.
        (see active_work.scde)

        Parameters
        ----------
        n : int
            Size of packs on which to average active work and order parameter.
        int_max : int or None
            Maximum number of packs consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of packs.
                  int_max cannot exceed the maximum number of nonoverlapping
                  packs.

        Returns
        -------
        hist : (*, 3) float Numpy array
            Values of the histogram:
                (0) Active work bin.
                (1) Order parameter bin.
                (2) Proportion.
            NOTE: This histogram is rather a probability density function,
                  therefore the integral over the bins is equal to 1 and thus
                  the values should be interpreted differently than a simple
                  proporition of observations.
        """

        workArray, orderParameter = self.nWorkOrder(n, int_max=int_max)
        pdf = PDF(workArray, orderParameter)

        return np.transpose(
            [*(lambda axes: [axes[:, -1], axes[:, -2]])(    # invert axes (work and order) order
                (pdf.extended_axes.reshape(np.prod(pdf.pdf.shape), 2))),
            pdf.pdf.flatten()])

    def meanStdCor(self, n=1, int_max=None):
        """
        Returns means anf standard deviations of active work and order
        parameter, and their Pearson correlation coefficient.

        Parameters
        ----------
        n : int
            Size of packs on which to average active work and order parameter.
        int_max : int or None
            Maximum number of packs consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of packs.
                  int_max cannot exceed the maximum number of nonoverlapping
                  packs.

        Returns
        -------
        meanWork : float
            Mean of active work.
        meanOrder : float
            Mean of order parameter.
        stdWork : float
            Standard deviation of active work.
        stdOrder : float
            Standard deviation of order parameter.
        corWorkOrder : float
            Pearson correlation coefficient of active work and order parameter.
        """

        workArray, orderParameter = self.nWorkOrder(n, int_max=int_max)

        meanWork = workArray.mean()
        meanOrder = orderParameter.mean()

        stdWork = workArray.std()
        stdOrder = orderParameter.std()

        corWorkOrder = np.cov(
            np.stack((workArray, orderParameter),
                axis=0))[0, 1]/(stdWork*stdOrder)

        return meanWork, meanOrder, stdWork, stdOrder, corWorkOrder
