"""
Module maths provides useful mathematic tools.

(taken from https://github.com/yketa/active_particles/tree/master/maths.py)
"""

import numpy as np

import math

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

def mean_sterr(values):
    """
    Returns mean and standard error of values.

    Parameters
    ----------
    values : float array
        Values.

    Returns
    -------
    mean : float
        Mean of values.
    sterr : float
        Standard error of values.
    """

    values = np.array(values)
    if values.size == 0: return None, None

    return np.mean(values), np.std(values)/np.sqrt(np.prod(values.shape))

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
