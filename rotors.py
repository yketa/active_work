"""
Module rotors provides classes and functions to compute and analyse
orientational dynamics and statistics of interacting Brownian rotors.

(see https://yketa.github.io/DAMTP_MSC_2019_Wiki/#N-interacting%20Brownian%20rotors)
"""

import numpy as np
from scipy.special import mathieu_a, mathieu_cem
import scipy.optimize as optimize

from active_work.read import DatR
from active_work.maths import Distribution

class Rotors(DatR):
    """
    Compute and analyse orientational dynamics and statistics from simulation
    data.

    (see https://yketa.github.io/DAMTP_MSC_2019_Wiki/#N-interacting%20Brownian%20rotors)
    """

    def __init__(self, filename, skip=1):
        """
        Loads file.

        Parameters
        ----------
        filename : string
            Name of input data file.
        skip : int
            Skip the `skip' first computed frames in the following calculations.
            (default: 1)
            NOTE: This can be changed at any time by setting self.skip.
        """

        super().__init__(filename)  # initialise with super class

        self.skip = skip    # skip the `skip' first frames in the analysis

    def nOrder(self, int_max=None, norm=False):
        """
        Returns array of order parameters.

        Parameters
        ----------
        int_max : int or None
            Maximum number of frames to consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of frames.
        norm : bool
            Return norm of order parameter rather than 2D order parameter.
            (default: False)

        Returns
        -------
        nu : [not(norm)] (*, self.N, 2) float numpy array
             [norm] (*, self.N) float numpy array
            Array of order parameters.
        """

        nu = []
        for time0 in self._time0(int_max=int_max):
            nu += [self.getOrderParameter(time0, norm=norm)]
        nu = np.array(nu)

        return nu

    def orderHist(self, nBins, int_max=None, vmin=0, vmax=1, log=False,
        rescaled_to_max=False):
        """
        Returns histogram with `nBins' bins of order parameter norm.

        Parameters
        ----------
        nBins : int
            Number of bins of the histogram.
        int_max : int or None
            Maximum number of frames to consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of frames.
        vmin : float
            Minimum value of the bins. (default: 0)
        vmax : float
            Maximum value of the bins. (default: 1)
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

        return Distribution(self.nOrder(int_max=int_max, norm=True)).hist(
                nBins, vmin=vmin, vmax=vmax, log=log,
                rescaled_to_max=rescaled_to_max)

    def nu_pdf_th(self, *nu):
        """
        Returns value of theoretical probability density function of the order
        parameter norm.

        (see https://yketa.github.io/DAMTP_MSC_2019_Wiki/#N-interacting%20Brownian%20rotors)

        Parameters
        ----------
        nu : float
            Values of the order parameter norm at which to evaluate the
            probability density function.

        Returns
        -------
        pdf : (*,) float numpy array
            Probability density function.
        """

        return nu_pdf_th(self.N, self.g, self.Dr, *nu)

    def _time0(self, int_max=None):
        """
        Returns array of frames at which to compute orientations.

        Parameters
        ----------
        int_max : int or None
            Maximum number of frames to consider. (default: None)
            NOTE: If int_max == None, then take the maximum number of frames.
                  WARNING: This can be very big.

        Returns
        -------
        time0 : (*,) float numpy array
            Array of frames.
        """

        if int_max == None: return np.array(range(self.skip, self.frames - 1))
        return np.linspace(
            self.skip, self.frames - 1, int(int_max), endpoint=False, dtype=int)

class Mathieu:
    """
    Provides estimates of the SCGF and the rate function of a single rotor from
    Mathieu functions, as well as optimal control potential for the angle.

    (see https://yketa.github.io/DAMTP_MSC_2019_Wiki/#Brownian%20rotors%20LDP)
    (see https://en.wikipedia.org/wiki/Mathieu_function)
    """

    def __init__(self, Dr):
        """
        Defines parameters.

        Parameters
        ----------
        Dr : float
            Rotational diffusivity.
        """

        self.Dr = Dr

        self._mathieu_order = 0 # order of Mathieu function

    def SCGF(self, *s):
        """
        Returns estimate of the SCGF.

        Parameters
        ----------
        s : float 2-uple
            Biasing parameter.

        Returns
        -------
        psi : float Numpy array
            Scaled cumulant generating function.
        """

        return np.array(list(map(
            lambda _s: -(self.Dr/4.)*(
                self._mathieu_characteristic_a(_s[0])
                + self._mathieu_characteristic_a(_s[1])),
            s)))

    def rate(self, *p):
        """
        Returns estimate of the rate function.

        Parameters
        ----------
        p : float 2-uple
            Polarisation norm.

        Returns
        -------
        I : float Numpy array
            Rate function.
        """

        return np.array(list(map(
            lambda _p: -np.array(optimize.minimize(
                lambda s: np.dot(s, _p) + self.SCGF(s)[0],
                (0, 0))['fun'], ndmin=1)[0],
            p)))

    def optimal_potential(self, s, *theta):
        """
        Returns optimal control potential for biasing parameter `s'.

        Parameters
        ----------
        s : float
            Biasing parameter.
        theta : float
            Angles (in radians) at which to evaluate the potential.

        Returns
        -------
        phi : float Numpy array
            Optimal control potential.
        """

        return (np.array(list(map(
            lambda _theta: -2*np.log(self._mathieu_function(s, _theta)),
            theta)))
            + 2*np.log(self._mathieu_function(s, 0)))   # normalisation

    def optimal_potential_curvature(self, s):
        """
        Returns curvature of optimal control potential at theta = 0.

        Parameters
        ----------
        s : float
            Biasing parameter.

        Returns
        -------
        k : float
            Curvature at theta = 0.
        """

        return (1./2.)*(
            self._mathieu_characteristic_a(s)
            - 2*self._mathieu_characteristic_q(s))

    def _mathieu_characteristic_q(self, s):
        """
        Returns characteristic value 'q' of the Mathieu function.

        Notation from https://en.wikipedia.org/wiki/Mathieu_function.

        Parameters
        ----------
        s : float
            Biasing parameter.

        Returns
        -------
        q : float
            Characteristic value 'q' of the Mathieu function.
        """

        return (2.*s)/self.Dr

    def _mathieu_characteristic_a(self, s):
        """
        Returns characteristic value 'a' of the Mathieu function.

        Notation from https://en.wikipedia.org/wiki/Mathieu_function.

        Parameters
        ----------
        s : float
            Biasing parameter.

        Returns
        -------
        a : float
            Characteristic value 'a' of the Mathieu function.
        """

        return mathieu_a(
            self._mathieu_order,
            self._mathieu_characteristic_q(s))

    def _mathieu_function(self, s, theta):
        """
        Returns Mathieu function evaluated at angle `theta

        Notation from https://en.wikipedia.org/wiki/Mathieu_function.

        Parameters
        ----------
        s : float
            Biasing parameter.
        theta : float
            Angle (in radians) at which to evaluate.

        Returns
        -------
        ce : float
            Mathieu function evaluated at `theta'.
        """

        return mathieu_cem(
            self._mathieu_order, self._mathieu_characteristic_q(s),
            (180./np.pi)*theta/2.)[0]

def nu_pdf_th(N, g, Dr, *nu):
    """
    Returns value of theoretical probability density function of the order
    parameter norm.

    (see https://yketa.github.io/DAMTP_MSC_2019_Wiki/#N-interacting%20Brownian%20rotors)

    Parameters
    ----------
    N : int
        Number of rotors.
    g : float
        Aligning torque parameter.
    Dr : float
        Rotational diffusivity.
    nu : float
        Values of the order parameter norm at which to evaluate the
        probability density function.

    Returns
    -------
    pdf : (*,) float numpy array
        Probability density function.
    """

    Z = (1 - np.exp(-N*(1 + g/Dr)))/(2*N*(1 + g/Dr))    # partition function

    return np.array(list(map(
        lambda _nu: _nu*np.exp(-N*(1 + g/Dr)*(_nu**2))/Z,
        nu)))
