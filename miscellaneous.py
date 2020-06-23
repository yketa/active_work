"""
Module miscellaneous contains miscellaneous functions and classes related to
this project.
"""

import numpy as np
import scipy.optimize as optimize

######################
### RTPs ON A RING ###
######################

class RTPring:
    """
    Theoretical results for 2 RTPs on a ring.

    (see https://github.com/yketa/DAMTP_MSC_2019_Wiki/tree/master/Summaries/RTPring)
    """

    def __init__(self):
        """
        Initialises computation parameters.
        """

        self._zero = 1e-16

    def Phi(self, Lambda):
        """
        Rescaled cumulant generating function (CGF).

        Parameters
        ----------
        Lambda : float
            Rescaled biasing parameter.

        Returns
        -------
        Phi : float
            Rescaled SCGF.
        """

        if Lambda == 0: return 0

        if Lambda > 0:
            return optimize.brentq(
                lambda _: self._scaling_function(_, Lambda),
                np.abs(self._zero),
                1e4)

        if Lambda < 0:
            return optimize.brentq(
                lambda _: self._scaling_function(_, Lambda),
                -((np.pi**2)/4 - np.abs(self._zero)),
                -np.abs(self._zero))

    def Gamma(self, Lambda):
        """
        Rescaled sticking probability of parallel particles.

        Parameters
        ----------
        Lambda : float
            Rescaled biasing parameter.

        Returns
        -------
        Gamma : float
            Rescaled sticking probability.
        """

        if Lambda == 0: return 1    # wild guess

        return self.Phi(Lambda)/Lambda

    def LEpsilon(self, Lambda, *rL):
        """
        Product of ring length and regular part of the probability density
        function.

        Parameters
        ----------
        Lambda : float
            Resclaed biasing parameter.
        rL : float
            Ratio of position over ring length.

        Returns
        -------
        Epsilon : float Numpy array
            Probability density function.
        """

        rL = np.array(rL)

        prefactor = (1./4.)*self.Gamma(Lambda)

        if Lambda == 0: return np.full(rL.shape, fill_value=prefactor)

        phi = self.Phi(Lambda)

        if Lambda > 0:
            return (prefactor*np.cosh(np.sqrt(phi)*(1 - 2*rL))
                /np.cosh(np.sqrt(phi)))
        if Lambda < 0:
            return (prefactor*np.cos(np.sqrt(-phi)*(1 - 2*rL))
                /np.cos(np.sqrt(-phi)))

    def _scaling_function(self, Phi, Lambda):
        """
        Function which root for a given `Lambda` gives the corresponding `Phi`.

        Parameters
        ----------
        Phi : float
            Rescaled cumulant generating function (CGF).
        Lambda : float
            Rescaled biasing parameter.

        Returns
        -------
        scale : float
            Scaling function.
        """

        if Lambda == 0: return 0

        if Lambda > 0:
            return (
                1./(np.tanh(np.sqrt(np.abs(Phi)))*np.sqrt(np.abs(Phi)))
                - 1./Lambda)
        if Lambda < 0:
            return (
                1./(np.tan(np.sqrt(np.abs(Phi)))*np.sqrt(np.abs(Phi)))
                + 1./Lambda)
