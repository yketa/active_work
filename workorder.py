"""
Module workorder provides classes to conjointly and analyse active work and
order parameter.
"""

import numpy as np

from active_work.work import ActiveWork
from active_work.maths import Histogram3D

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

    # def corWorkOrderIns(self,
    #     tau0=1, n_max=100, int_max=None, min=None, max=None, log=True):
    #     """
    #     Compute correlations of the fluctuations of the active work and the
    #     order parameter averaged over `tau0' between different times.
    #
    #     Parameters
    #     ----------
    #     tau0 : int
    #         Number of consecutive individual active works on which to average
    #         it. (default: 1)
    #     n_max : int
    #         Maximum number of values at which to evaluate the correlation.
    #         (default: 100)
    #     int_max : int or None
    #         Maximum number of different intervals to consider in order to
    #         compute the mean which appears in the correlation expression.
    #         (default: None)
    #         NOTE: if int_max == None, then a maximum number of disjoint
    #               intervals will be considered.
    #     min : int or None
    #         Minimum value at which to compute the correlation. (default: None)
    #         NOTE: if min == None then min = `tau0', otherwise the minimum of
    #               `tau0' and `min' is taken.
    #     max : int or None
    #         Maximum value at which to compute the correlation. (default: None)
    #         NOTE: this value is passed as `max' to self.n.
    #     log : bool
    #         Logarithmically space values at which the correlations are
    #         computed. (default: True)
    #
    #     Returns
    #     -------
    #     cor : (5, *) numpy array
    #         Array of:
    #             (0) value at which the correlation is computed,
    #             (1) mean of the computed correlation,
    #             (2) standard error of the computed correlation,
    #             (3) standard deviation of the active work computed at the
    #                 beginning of the interval,
    #             (4) standard deviation of the active work computed at the end
    #                 of the interval.
    #     """

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
