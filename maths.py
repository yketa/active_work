"""
Module maths provides useful mathematic tools.

(see also https://github.com/yketa/active_particles/tree/master/maths.py)
"""

import numpy as np
from numpy.polynomial.polynomial import polyadd, polypow
import math
from collections import OrderedDict
from scipy import optimize

from active_work.scde import PDF

#####################
### MISCELLANEOUS ###
#####################

def relative_positions(positions, point, box_size):
    """
    Returns relative positions to point in box of extent
    (-box_size/2, box_size) in both dimensions of space.

    Parameters
    ----------
    positions : float array
        Position of single point or array of positions.
    point : float array
        Position of the new centre.
    box_size : float or array
        Length of the box in one dimension or all dimensions.

    Returns
    -------
    rel_positions : float array
        Relative positions.
    """

    return (np.array(positions) - np.array(point)
        + np.array(box_size)/2)%np.array(box_size) - np.array(box_size)/2

def wo_mean(arr):
    """
    Returns deviation of values in array with respect to mean of array.

    Parameters
    ----------
    arr : array like
        Array of values.

    Returns
    -------
    dev_arr : array like
        Deviations from mean of array.
    """

    return np.array(arr) - np.mean(arr, axis=0)

class DictList(dict):
    """
    Custom hash table class to give value [] to uninitialised keys.
    """
    def __init__(self):
        super().__init__()
    def __getitem__(self, key):
        try:
            return super().__getitem__(key)
        except KeyError:
            return []

def normalise1D(*vector):
    """
    Returs 1D vector of unitary norm with same direction.

    Parameters
    ----------
    vector : 1D array-like or coordinates as positional arguments
        Vector to normalise.

    Returns
    -------
    u_vector : 1D Numpy array
        Unitary vector with same direction.
    """

    vector = np.array(vector).flatten() # 1D vector

    norm = np.linalg.norm(vector)   # vector norm
    if norm == 0: return vector     # vector is 0
    return vector/norm

def amplogwidth(arr, factor=2):
    """
    Calculates the amplitudes of elements in array arr and, excluding the
    zeros, returns the mean of the logarithms of these amplitudes plus and
    minus factor times their standard deviation.

    Parameters
    ----------
    arr : array like
        Array.
    factor : float
        Width factor. (default: 2)

    Returns
    -------
    min : float
        E(log(||arr||)) - factor*V(log(||arr||))
    max : float
        E(log(||arr||)) + factor*V(log(||arr||))
    """

    log = np.ma.log10(np.sqrt(np.sum(arr**2, axis=-1))) # logarithms of amplitudes
    mean = log.mean()                                   # means of logarithms of amplitudes
    std = log.std()                                     # standard deviation of logarithms of amplitudes

    return mean - factor*std, mean + factor*std

def mean_sterr(values, remove=False, max=None):
    """
    Returns mean and standard error of values.

    Parameters
    ----------
    values : float array
        Values.
    remove : bool
        Remove inf and -inf as well as nan. (default: False)
        NOTE: A warning will be issued if remove == False and such objects are
              encountered.
        NOTE: This is not guaranteed to work non-1D arrays as the shape may
              change.
    max : float or None
        Remove data which is strictly above max in absolute value.
        (default: None)
        NOTE: max != None will trigger remove = True.

    Returns
    -------
    mean : float or float Numpy array
        Mean of values.
    sterr : float or float Numpy array
        Standard error of values = std(...)/sqrt(len(...)).
        NOTE: This is relevant if all values are independent.
    """

    if max != None: remove = True

    values = np.array(values)
    if remove: values = (
        (lambda _: _[~np.isinf(_)])(    # remove inf
        (lambda __: __[~np.isnan(__)])( # remove nan
        values)))
    if max != None: values = values[np.abs(values) <= max]
    if values.size == 0: return None, None

    return values.mean(axis=0), values.std(axis=0)/np.sqrt(values.shape[0])

def aggregate(array, n):
    """
    Aggregate `array' by taking mean of bits of `n' consecutive values on last
    axis.

    Parameters
    ----------
    array : float array-like
        Array to aggregate.
    n : int
        Size of bits to aggregate.

    Returns
    -------
    out : float Numpy array
        Aggregated array.
    """

    array = np.array(array)

    return (np.array(
        np.split(array, array.shape[-1]//n, axis=-1))
        .transpose(
            tuple(i for i in range(1, len(array.shape)))
            + (0, len(array.shape)))
        .mean(axis=-1))

def cov(array1, array2):
    """
    Return covariance of variables `array1' and `array2'.

        Cov(`array1', `array2') = <`array1'*`array2'> - <`array1'><`array2'>

    NOTE: `array1' and `array2' have to be the same size.

    Parameters
    ----------
    array1 : float array-like
        First set of measurements.
    array2 : float array-like
        Second set of measurements.

    Returns
    -------
    covariance : float
        Cov(`array1', `array2') = <`array1'*`array2'> - <`array1'><`array2'>.
    """

    return (lambda a, b: (a*b).mean() - a.mean()*b.mean())(
        *(np.array(array1), np.array(array2)))

def divide_arrays(array1, array2):
    """
    Divide array1 by array2, and outputs 0 values where array2 is equal to 0.
    NOTE: array1, array2 and out must have the same shapes.

    Parameters
    ----------
    array1 : array-like
        Numerator array.
    array2 : array-like
        Denominator array.

    Returns
    -------
    array : array-like
        Quotient array.
    """

    if not(isinstance(array1, np.ndarray)): array1 = np.array(array1)
    if not(isinstance(array2, np.ndarray)): array2 = np.array(array2)

    return np.divide(array1, array2,
        out=np.zeros(array1.shape, dtype=array1.dtype), where=array2!=0)

def linspace(init, fin, number, endpoint=True):
    """
    Returns linearly spaced integer between `init' and `fin' with a maximum of
    `number' of them.

    Parameters
    ----------
    init : int
        Minimum value.
    fin : int
        Maximum value.
    number : int
        Number of values.
    endpoint : bool
        Include `number' in the array.

    Returns
    -------
    values : numpy array
        Array of values.
    """

    return np.array(list(OrderedDict.fromkeys(np.linspace(
        init, fin, number,
        endpoint=endpoint, dtype=int))))

def logspace(init, fin, number, endpoint=True):
    """
    Returns logarithmically spaced integer between `init' and `fin' with a
    maximum of `number' of them.

    Parameters
    ----------
    init : int
        Minimum value.
    fin : int
        Maximum value.
    number : int
        Number of values.
    endpoint : bool
        Include `number' in the array.

    Returns
    -------
    values : numpy array
        Array of values.
    """

    return np.array(list(OrderedDict.fromkeys(map(lambda x: int(round(x)),
        np.exp(np.linspace(
            np.log(init), np.log(fin), number,
            endpoint=endpoint))))))

def meanStdCut(array, cut=None):
    """
    Returns mean and standard deviation of array with values farther than
    `cut' * array.std() from the mean removed.

    Parameters
    ----------
    array : array-like
        Array of values.
    cut : float
        Width in units of array.std() to consider. (default: None)
        NOTE: if cut==None, then no value is excluded.

    Returns
    -------
    mean : float
        Mean of the truncated ensemble.
    std : float
        Standard deviation of the truncated ensemble.
    """

    array = np.array(array)

    if cut == None: return array.mean(), array.std()

    array = array[np.abs(array - array.mean()) < cut*array.std()]
    return array.mean(), array.std()

def angle(dx, dy):
    """
    Returns angle from in x- and y-coordinates.

    Parameters
    ----------
    dx : float
        x-coordinate or difference in x-coordinate.
    dy : float
        y-coordinate or difference in y-coordinate.

    Returns
    -------
    ang : float
        Corresponding angle in radians.
    """

    return math.atan2(dy, dx)

class CurveFit:
    """
    Fit a 1-variable scalar function to data points and evaluate with
    uncertainty. (see scipy.optimize.curve_fit)
    """

    def __init__(self, func, xdata, ydata, jac=None, **kwargs):
        """
        Fit curve to data points. (see scipy.optimize.curve_fit)

        Parameters
        ----------
        func : callable scalar function
            Model function, f(x, ...), taking the independent variable as the
            first argument and the parameters to fit as separate remaining
            arguments.
        xdata : float array-like
            x-data to fit.
        ydata : float array-like
            y-data to fit.
        jac : callable 1D-array-like function or None
            Jacobian matrix of the model function with respect to parameters.
            (default: None)
            NOTE: standard deviations will not be computed if jac == None.

        Optional keyword arguments
        --------------------------
        Additional keyword arguments will be passed to scipy.optimize.curve_fit.
        """

        self.func = func
        self.jac = jac
        self.xdata = np.array(xdata)
        self.ydata = np.array(ydata)
        self.curve_fit_kwargs = kwargs

        self.popt, self.pcov = optimize.curve_fit(
            self.func, self.xdata, self.ydata,
            #jac=self.jac,
            **self.curve_fit_kwargs)

    def eval(self, *x):
        """
        Evaluate fitted curve with standard deviation.

        Parameters
        ----------
        x : float
            x-data.

        Returns
        -------
        y : float Numpy array
            Fitted y-data.
        sigma : float Numpy array
            Uncertainty on y-data.
        """

        y = np.array(list(map(lambda _: self.func(_, *self.popt), x)))
        sigma = np.array(list(map(lambda _: self._sigma(_), x)))

        return y, sigma

    def _sigma(self, x):
        """
        Evaluate standard deviation at point `x'.

        Parameters
        ----------
        x : float
            x-data.

        Returns
        -------
        sigma : float
            Standard deviation at `x'.
            NOTE: returns 0 if self.jac == None.
        """

        if self.jac is None: return 0

        return np.sqrt(
            np.dot(
                self.jac(x, *self.popt),
                np.dot(
                    self.pcov,
                    np.transpose(self.jac(x, *self.popt)))))

###################
### POLYNOMIALS ###
###################

def evalPol(pol, *x):
    """
    Evaluate polynomial. (see numpy.polyeval)

    NOTE: Returns None if x == None.

    Parameters
    ----------
    pol : (*,) float array-like
        Polynomial coefficients (highest first).
    x : float
        Values at which to evaluate the polynomial.

    Returns
    -------
    y : (*,) float Numpy array
        Evaluated polynomial.
    """

    return np.array(list(map(
        lambda _x: None if _x == None else np.polyval(pol, _x),
        x)))

def addPol(*pol):
    """
    Wrapper around numpy.polynomial.polynomial.polyadd to sum several
    polynomials at once.

    Parameters
    ----------
    pol : (*,) float array-like
        Polynomial coefficients (lowest first).

    Returns
    -------
    sum : (*,) float array-like
        Sum of polynomials.
    """

    sum = np.zeros((1,))
    for p in pol:
        sum = polyadd(sum, p)

    return sum

class Polynomial:
    """
    Store and evaluate polynomials with covariance matrix on their
    coefficients.
    """

    def __init__(self, pol, cov):
        """
        Compute standard devation polynomial.

        Considering we have polynomial

            P(x) = \\sum_i a_i x^i,

        and \\sigma_ij the covariance of the i-th and j-th coefficients, a_i and
        a_j, we compute

            \\sigma_f(x) = \\sqrt(\\sum_ij \\sigma_ij x^{i+j})

        as the variance of the polynomial f evaluated at x.

        (see https://en.wikipedia.org/wiki/Propagation_of_uncertainty)

        Parameters
        ----------
        pol : (N,) float array-like
            Polynomial coefficients (highest first).
        cov : (N, N) float-array like
            Covariance matrix.
        """

        pol, cov = np.array(pol), np.array(cov)

        if cov.shape != (pol.size, pol.size): raise ValueError

        self.pol = pol
        self.deg = len(self.pol) - 1

        self.cov = cov
        self.covpol = np.zeros((2*self.deg + 1,))   # covariance polynomial for estimation of standard deviation
        for i in range(self.deg + 1):
            for j in range(self.deg + 1):
                self.covpol[i + j] += self.cov[i, j]

    def eval(self, *x, sigma=False):
        """
        Evaluate polynomial and standard deviation.

        Parameters
        ----------
        x : float
            x-values at which to evaluate the polynomial.
        sigma : bool
            Return associated standard deviation. (default: False)

        Returns
        -------
        y : (*,) float Numpy array
            Evaluated polynomial.
        std : [if sigma] (*,) float Numpy array
            Evaluated standard deviation.
        """

        y = evalPol(self.pol, *x)
        if sigma: return y, np.sqrt(evalPol(self.covpol, *x))
        return y

class CompPol(Polynomial):
    """
    Create polynomial from composition of two others.
    """

    def __init__(self, pol1, pol2):
        """
        Composes two polynomials (i.e., `pol1'(`pol2')) with distinct covariance
        matrices.

        Considering we have polynomials

            f(x) = \\sum_i a_i x^i,
            g(x) = \\sum_j b_j x^j,

        with variances \\sigma_f and \\sigma_g when evaluated (see
        active_work.maths.Polynomial), we compute

            \\sigma_fg(x) = \\sigma_f(g(x))
                + \\sigma_g(x) \\times [ \\sum_i i a_i g(x)^{i-1} ]^2,

        as the variance of the composed polynomial f(g) evaluated at x. We
        stress that this considers no correlations between the coefficients of
        the polynomials.

        (see https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Non-linear_combinations)

        Parameters
        ----------
        pol1 : active_work.maths.Polynomial
            First polynomial.
        pol2 : active_work.maths.Polynomial
            Second polynomial.

        Returns
        -------
        pol : active_work.maths.Polynomial
            Composed polynomial.
        """

        self._pol1, self._pol2 = pol1, pol2
        self.deg = self._pol1.deg*self._pol2.deg    # degree of composed polynomial

        # WARNING: numpy.polynomial.polynomial.polyadd and polypow considers
        #          arrays as polynomials with lowest coefficient first,
        #          contrarily to polyval and polyfit.
        _pol1, _pol2 = self._pol1.pol[::-1], self._pol2.pol[::-1]

        self.pol = np.zeros((1,))   # composed polynomial
        for i in range(pol1.deg + 1):
            self.pol = polyadd(self.pol, _pol1[i]*polypow(_pol2, i))

        self.pol = self.pol[::-1]

    def eval(self, *x, sigma=False):
        """
        Evaluate polynomial and standard deviation.

        Parameters
        ----------
        x : float
            x-values at which to evaluate the polynomial.
        sigma : bool
            Return associated standard deviation. (default: False)

        Returns
        -------
        y : (*,) float Numpy array
            Evaluated polynomial.
        std : [if sigma] (*,) float Numpy array
            Evaluated standard deviation.
        """

        y = self._pol1.eval(*self._pol2.eval(*x, sigma=False), sigma=False)

        if sigma:
            var = np.array(list(map(
                lambda _x, _gx:
                    evalPol(self._pol1.covpol, _gx)[0]
                    + evalPol(self._pol2.covpol, _x)[0]*(np.sum(
                        [k*self._pol1.pol[-1-k]*(_gx**(k - 1))
                            for k in range(1, self._pol1.deg + 1)])**2),
                *(x, self._pol2.eval(*x, sigma=False)))))
            return y, np.sqrt(var)

        return y

class PolyFit(Polynomial):
    """
    Perform least squares polynomial fit and evaluate fit. (see numpy.polyfit)
    """

    def __init__(self, x, y, deg=1):
        """
        Perform fit.

        Parameters
        ----------
        x : (*,) float array-like
            x-coordinates of sample points.
        y : (*,) float  array-like
            y-coordinates of sample points.
        deg : int
            Degree of the fitting polynomial. (default: 1)
        """

        self.x = np.array(x)
        self.xmin = self.x.min()
        self.xmax = self.x.max()

        self.y = np.array(y)
        self.ymin = self.y.min()
        self.ymax = self.y.max()

        pol, cov = np.polyfit(self.x, self.y, deg, cov=True)
        Polynomial.__init__(self, pol, cov)

    def eval(self, *x, restrict=False, sigma=False):
        """
        Evaluate polynomial and standard deviation.
        (see active_work.maths.Polynomial.eval)

        Parameters
        ----------
        x : float
            x-values at which to evaluate the polynomial.
        restrict : bool
            Discard x-values which are not in the range of the original data.
            (default: False)
        sigma : bool
            Return associated standard deviation. (default: False)

        Returns
        -------
        y : (*,) float Numpy array
            Evaluated polynomial.
        std : [if sigma] (*,) float Numpy array
            Evaluated standard deviation.
        """

        if restrict:
            x = np.array(x)
            for i in range(len(x)):
                if x[i] == None: continue
                if x[i] < self.xmin or x[i] > self.xmax: x[i] = None

        return super().eval(*x, sigma=sigma)

#####################
### DISTRIBUTIONS ###
#####################

class Distribution:
    """
    Analyse distribution from array of values.
    """

    def __init__(self, valuesArray):
        """
        Define array of values.

        Parameters
        ----------
        valuesArray : float array-like
            Array of values.
        """

        self.valuesArray = np.array(valuesArray).flatten()

        self.min = self.valuesArray.min()
        self.max = self.valuesArray.max()

        self.mean = self.valuesArray.mean()
        self.std = self.valuesArray.std()

    def pdf(self):
        """
        Returns probability density function from array of values.

        Returns
        -------
        axes : numpy array
            Values at which the probability density function is evaluated.
        pdf : float numpy array
            Values of the probability density function.
        """

        pdf = PDF(self.valuesArray)
        return pdf.axes[0], pdf.pdf

    def hist(self, nBins, vmin=None, vmax=None, log=False,
        rescaled_to_max=False):
        """
        Returns histogram of array of values.

        Parameters
        ----------
        nBins : int
            Number of bins of the histogram.
        vmin : float
            Minimum value of the bins. (default: None)
            NOTE: if vmin==None, then minimum of array is taken.
        vmax : float
            Maximum value of the bins. (default: None)
            NOTE: if vmax==None, then maximum of array is taken.
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

        if vmin == None: vmin = self.min
        if vmax == None: vmax = self.max
        histogram = Histogram(nBins, vmin, vmax, log=False)
        histogram.values = self.valuesArray

        bins = histogram.bins
        hist = histogram.get_histogram()
        if rescaled_to_max: hist /= hist.max()
        if not(log): return bins, hist
        else: return bins[hist > 0], np.log(hist[hist > 0])

    def gauss(self, *x, cut=None, rescaled_to_max=False):
        """
        Returns values of the Gaussian function corresponding to the mean and
        variance of the array of values.

        Parameters
        ----------
        x : float
            Values at which to evaluate the Gaussian function.
        cut : float or None
            Width in units of the standard deviation of the array of values to
            consider when computing mean and standard deviation.
            (see self._meanStdCut) (default: None)
            NOTE: if cut==None, the width is taken to infinity, i.e. no value is
                  excluded.
        rescaled_to_max : bool
            Rescale function by its computed maximum. (default: False)

        Returns
        -------
        gauss : float numpy array
            Values of the Gaussian function at x.
        """

        mean, std = self._meanStdCut(cut=cut)

        if rescaled_to_max: norm = 1
        else: norm = np.sqrt(2*np.pi*(std**2))

        gauss = lambda y: (
            np.exp(-((y - mean)**2)/(2*(std**2)))
            /norm)

        return np.array(list(map(gauss, x)))

    def _meanStdCut(self, cut=None):
        """
        Returns mean and standard deviation of values of array with values
        farther than `cut' * self.valuesArray.std() if the mean removed.

        Parameters
        ----------
        array : array-like
            Array of values.
        cut : float
            Width in units of self.valuesArray.std() to consider.
            (default: None)
            NOTE: if cut==None, then no value is excluded.

        Returns
        -------
        mean : float
            Mean of the truncated ensemble.
        std : float
            Standard deviation of the truncated ensemble.
        """

        return meanStdCut(self.valuesArray, cut=cut)

class JointDistribution:
    """
    Analyse joint distribution from 2 arrays of values.
    """

    def __init__(self, valuesArray1, valuesArray2):
        """
        Define array of values.

        Parameters
        ----------
        valuesArray1 : float array-like
            First array of values.
        valuesArray2 : float array-like
            Second array of values.
        """

        self.valuesArray1 = np.array(valuesArray1).flatten()
        self.valuesArray2 = np.array(valuesArray2).flatten()

        self.min1 = self.valuesArray1.min()
        self.max1 = self.valuesArray1.max()
        self.min2 = self.valuesArray2.min()
        self.max2 = self.valuesArray2.max()

    def pdf(self):
        """
        Returns joint probability density function from arrays of values.

        Returns
        -------
        pdf3D : (*, 3) float Numpy array
            (0) Value of the first quantity at which the PDF is evaluated.
            (1) Value of the second quantity at which the PDF is evaluated.
            (2) PDF.
        """

        pdf = PDF(self.valuesArray1, self.valuesArray2)

        return np.transpose(
            [*(lambda axes: [axes[:, -1], axes[:, -2]])(    # invert axes order
                (pdf.extended_axes.reshape(np.prod(pdf.pdf.shape), 2))),
            pdf.pdf.flatten()])

    def hist(self, nBins, vmin1=None, vmax1=None, vmin2=None, vmax2=None):
        """
        Returns 3D histogram of arrays of values.

        Parameters
        ----------
        nBins : int or 2-uple-like of int
            Number of bins of the histogram in all or each direction.
        vmin1 : float
            Minimum value of the bins for the first array. (default: None)
            NOTE: if vmin1==None, then minimum of array is taken.
        vmax1 : float
            Maximum value of the bins for the first array. (default: None)
            NOTE: if vmax1==None, then maximum of array is taken.
        vmin2 : float
            Minimum value of the bins for the second array. (default: None)
            NOTE: if vmin2==None, then minimum of array is taken.
        vmax2 : float
            Maximum value of the bins for the second array. (default: None)
            NOTE: if vmax2==None, then maximum of array is taken.

        Returns
        -------
        hist : (nBins.prod(), 3) float Numpy array
            Values of the histogram:
                (0) Bin value of the first quantity.
                (1) Bin value of the second quantity.
                (2) Proportion.
        """

        if vmin1 == None: vmin1 = self.min1
        if vmax1 == None: vmax1 = self.max1
        if vmin2 == None: vmin2 = self.min2
        if vmax2 == None: vmax2 = self.max2
        histogram = Histogram3D(nBins, (vmin1, vmin2), (vmax1, vmax2),
            log=False)
        histogram.values = np.array(list(
            zip(self.valuesArray1, self.valuesArray2)))

        return histogram.get_histogram()

##################
### HISTOGRAMS ###
##################

class Histogram:
    """
    Make histogram from lists of float values.
    """

    def __init__(self, Nbins, vmin, vmax, log=False):
        """
        Parameters
        ----------
        Nbins : int
            Number of histogram bins.
        vmin : float
            Minimum included value for histogram bins.
            NOTE: values lesser than vmin will be ignored.
        vmax : float
            Maximum excluded value for histogram bins.
            NOTE: values greater or equal to vmax will be ignored.
        log : bool.
            Logarithmically spaced histogram values. (default: False)
        """

        self.Nbins = int(Nbins)
        self.vmin = vmin
        self.vmax = vmax

        if log:
            self.bins = np.logspace(np.log10(self.vmin), np.log10(self.vmax),
                self.Nbins, endpoint=False, base=10)    # histogram bins
        else:
            self.bins = np.linspace(self.vmin, self.vmax,
                self.Nbins, endpoint=False)             # histogram bins

        self.reset_values()                 # reset values from which to compute the histogram
        self.hist = np.empty(self.Nbins)    # values of the histogram at bins

    def add_values(self, *values, replace=False):
        """
        Add values from which to compute the histogram.

        Parameters
        ----------
        values : float or float array-like
            Values to add.
        replace : bool
            Replace existing values. (default: False)
        """

        if replace: self.reset_values()
        for value in values: self.values = np.append(self.values, value)

    def reset_values(self):
        """
        Delete values from which to compute the histogram (self.values).
        """

        self.values = np.array([])

    def get_histogram(self):
        """
        Get histogram from values in self.values.

        Returns
        -------
        hist : Numpy array
            Values of the histogram at self.bins.
        """

        for bin in range(self.bins.size):
            bin_inf = self.bins[bin]
            try: bin_sup = self.bins[bin + 1]
            except IndexError: bin_sup = self.vmax
            self.hist[bin] = np.sum(
                (self.values >= bin_inf)*(self.values < bin_sup))

        binned_values = np.sum(self.hist)
        if binned_values == 0: return self.hist # no binned value
        else: self.hist /= np.sum(self.hist)
        return self.hist

class Histogram3D:
    """
    Make 3D histogram from lists of float 2-uples-like.
    """

    def __init__(self, Nbins, vmin, vmax, log=False):
        """
        Parameters
        ----------
        Nbins : int or 2-uple-like of int
            Number of histogram bins in each direction.
        vmin : float or 2-uple like of float
            Minimum included value for histogram bins.
            NOTE: values lesser than vmin will be ignored.
        vmax : float or 2-uple like of float
            Maximum excluded value for histogram bins.
            NOTE: values greater or equal to vmax will be ignored.
        log : bool.
            Logarithmically spaced histogram values. (default: False)
        """

        Nbins = np.array(Nbins, ndmin=1, dtype=int)
        self.Nbins = np.array([Nbins[0], Nbins[-1]])

        vmin, vmax = np.array(vmin, ndmin=1), np.array(vmax, ndmin=1)
        self.vmin = np.array([vmin[0], vmin[-1]])
        self.vmax = np.array([vmax[0], vmax[-1]])

        self.bins = []
        for dim in range(2):
            if log:
                self.bins += [np.logspace(
                    np.log10(self.vmin[dim]), np.log10(self.vmax[dim]),
                    self.Nbins[dim], endpoint=False, base=10)]  # histogram bins
            else:
                self.bins += [np.linspace(
                    self.vmin[dim], self.vmax[dim],
                    self.Nbins[dim], endpoint=False)]           # histogram bins

        self.reset_values()                             # reset values from which to compute the histogram
        self.hist = np.empty((self.Nbins.prod(), 3))    # values of the histogram at bins
        for bin0 in range(self.bins[0].size):
            self.hist[
                bin0*self.bins[1].size:(bin0 + 1)*self.bins[1].size, 0] = (
                self.bins[0][bin0])
            self.hist[
                bin0*self.bins[1].size:(bin0 + 1)*self.bins[1].size, 1] = (
                self.bins[1])

    def add_values(self, *values, replace=False):
        """
        Add values from which to compute the histogram.

        Parameters
        ----------
        values : float or float array-like
            Values to add.
        replace : bool
            Replace existing values. (default: False)
        """

        if replace: self.reset_values()
        for value in values: self.values += [tuple(value)]

    def reset_values(self):
        """
        Delete values from which to compute the histogram (self.values).
        """

        self.values = []

    def get_histogram(self):
        """
        Get histogram from values in self.values.

        Returns
        -------
        hist : (self.Nbins.prod(), 3) float Numpy array
            Values of the histogram:
                (0) Value of first axis bin.
                (1) Value of second axis bin.
                (2) Proportion.
        """

        values_array = np.array(self.values)
        for bin0 in range(self.bins[0].size):
            bin_inf0 = self.bins[0][bin0]
            try: bin_sup0 = self.bins[0][bin0 + 1]
            except IndexError: bin_sup0 = self.vmax[0]
            values = values_array[
                (values_array[:, 0] >= bin_inf0)
                *(values_array[:, 0] < bin_sup0)][:, 1]
            for bin1 in range(self.bins[1].size):
                bin_inf1 = self.bins[1][bin1]
                try: bin_sup1 = self.bins[1][bin1 + 1]
                except IndexError: bin_sup1 = self.vmax[1]
                self.hist[bin0*self.Nbins[1] + bin1, 2] = (
                    np.sum((values >= bin_inf1)*(values < bin_sup1)))

        if np.sum(self.hist[:, 2]) > 0: # there are binned values
            self.hist[:, 2] /= np.sum(self.hist[:, 2])
        return self.hist

#############
### GRIDS ###
#############

def vector_vector_grid(vector1, vector2, dtype=None):
    """
    From vector1 = (v1_i)_i and vector2 = (v2_i)_i, returns matrix
    M = (M_{i, j})_{i, j} = ((v1_i, v2_j))_{i, j}.

    Parameters
    ----------
    vector1 : 1D array-like
        Vector 1.
    vector2 : 1D array-like
        Vector 2.
    dtype : Numpy array dtype
        Data type of the Numpy array to return. (default: None)
        NOTE: if dtype == None, then the array is not converted to any type.

    Returns
    -------
    M : 2D array-like
        Matrix M.
    """

    M = np.zeros((len(vector1), len(vector2), 2))
    M[:, :, 0] = vector1
    M = np.transpose(M, (1, 0, 2))
    M[:, :, 1] = vector2

    if dtype != None: return M.astype(dtype)
    else: return M

def wave_vectors_2D(nx, ny, d=1):
    """
    Returns wave vectors for 2D signals with window lengths nx and ny in the
    two directions and sample spacing d.

    Parameters
    ----------
    nx : int
        Window length in first direction.
    ny : int
        Window length in second direction.
    d : float
        Sample spacing. (default: 1)

    Returns
    -------
    wave_vectors : (nx, ny, 2) Numpy array
        Grid of wave vectors.
    """

    return 2*np.pi*vector_vector_grid(
        np.fft.fftfreq(nx, d=d),
        np.fft.fftfreq(ny, d=d))

def g2Dto1D(g2D, L):
    """
    Returns cylindrical average of 2D grid.

    Parameters
    ----------
    g2D : 2D array
        2D grid.
        NOTE: g2D[0, 0] is considered the r=0 point on the grid, and we
        consider periodic boundaries.
    L : float or float array
        Length of the box represented by the grid in one dimension or all
        dimensions.

    Returns
    -------
    g1D : Numpy array
        Array of (r, g1D(r)) with g1D(r) the averaged 2D grid at radius r.
    """

    g2D = np.array(g2D)
    dL = np.array(L)/np.array(g2D.shape)    # boxes separation in each direction
    r_max = np.min(L)/2                     # maximum radius to be calculated in number of boxes

    g1D_dic = DictList()    # hash table of radii and values at radii

    for i in range(g2D.shape[0]):
        for j in range(g2D.shape[1]):
            radius = np.sqrt(np.sum((np.array((i, j))*dL)**2))  # radius corresponding to coordinates [i, j], [-i, j], [i, -j], [-i, -j]
            if radius <= r_max:
                g1D_dic[radius] += [g2D[i, j], g2D[-i, j], g2D[i, -j],
                    g2D[-i, -j]]

    return np.array(list(map(
        lambda radius: [radius, np.mean(g1D_dic[radius])],
        sorted(g1D_dic))))

def g2Dto1Dgrid(g2D, grid, average_grid=False):
    """
    Returns cylindrical average of square 2D grid with values of radius given
    by other parameter grid.

    Parameters
    ----------
    g2D : 2D array
        Square 2D grid.
    grid : 2D array
        Array of radii.
    average_grid : bool
        Return g2D grid with cylindrically averaged values.

    Returns
    -------
    g1D : Numpy array
        Array of (r, g1D(r)) with g1D(r) the averaged 2D grid at radius r.
    g2D_cylindrical [average_grid] : Numpy array
        Cylindrically averaged g2D.
    """

    g2D = np.array(g2D)
    grid = np.array(grid)

    g1D_dic = DictList()    # hash table of radii and values at radii

    for i in range(g2D.shape[0]):
        for j in range(g2D.shape[1]):
            g1D_dic[grid[i, j]] += [g2D[i, j]]

    g1D = np.array(list(map(
        lambda radius: [radius, np.mean(g1D_dic[radius])],
        sorted(g1D_dic))))

    if not(average_grid): return g1D

    g2D_cylindrical = np.zeros(grid.shape)
    for radius, mean_g in zip(*np.transpose(g1D)):
        for i, j in zip(*np.where(grid == radius)):
            g2D_cylindrical[i, j] = mean_g

    return g1D, g2D_cylindrical
