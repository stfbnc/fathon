#    mfdcca.pyx - mfdcca algorithm of fathon package
#    Copyright (C) 2019-  Stefano Bianchi
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
import numpy as np
from cython.parallel import prange
import ctypes
import pickle
import warnings
class MFDCCA:
    """MultiFractal Detrended Cross-Correlation Analysis class.

    Parameters
    ----------
    tsVec1 : iterable
        First time series used for the analysis.
    tsVec2 : iterable
        Second time series used for the analysis.
    """

    def __init__(self, tsVec1, tsVec2=[]):
        pass

    def computeFlucVec(self, winSizes, qList, polOrd=1, revSeg=False):
        """Computation of the fluctuations in each window for each q-order.

        Parameters
        ----------
        winSizes : numpy ndarray
            Array of window's sizes.
        qList : float or iterable or numpy ndarray
            List of q-orders used to compute `F`.
        polOrd : int, optional
            Order of the polynomial to be fitted in each window (default : 1).
        revSeg : bool, optional
            If True, the computation of `F` is repeated starting from the end
            of the time series (default : False).

        Returns
        -------
        numpy ndarray
            Array `n` of window's sizes.
        numpy ndarray
            qxn array `F` containing the values of the fluctuations in each
            window for each q-order.
        """
        return 0

    def fitFlucVec(self, nStart=-999, nEnd=-999, logBase=np.e, verbose=False):
        """Fit of the fluctuations values.

        Parameters
        ----------
        nStart : int, optional
            Size of the smaller window used to fit `F` at each q-order (default : first value of `n`).
        nEnd : int, optional
            Size of the bigger window used to fit `F` at each q-order (default : last value of `n`).
        logBase : float, optional
            Base of the logarithm for the log-log fit of `n` vs `F` (default : e).
        verbose : bool, optional
            Verbosity (default : False).

        Returns
        -------
        numpy ndarray
            Slope of the fit for each q-order.
        numpy ndarray
            Intercept of the fit for each q-order.
        """
        return 0

    def computeMassExponents(self):
        """Computation of the mass exponents.

        Returns
        -------
        numpy ndarray
            Mass exponents.
        """
        return 0

    def computeMultifractalSpectrum(self):
        """Computation of the multifractal spectrum.

        Returns
        -------
        numpy ndarray
            Singularity strengths.
        numpy ndarray
            Multifractal spectrum.
        """
        return 0

    def saveObject(self, outFileName):
        """Save current object state to binary file.
        
        Parameters
        ----------
        outFileName : str
            Output binary file. `.fathon` extension will be appended to the file name.
        """
        return 0

