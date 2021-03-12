#    dcca.pyx - dcca algorithm of fathon package
#    Copyright (C) 2019-2021  Stefano Bianchi
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
from . import dfa
class DCCA:
    """Detrended Cross-Correlation Analysis class.

    Parameters
    ----------
    n : numpy ndarray
        Array of window's sizes used for the computation.
    tsVec1 : iterable
        First time series used for the analysis.
    tsVec2 : iterable
        Second time series used for the analysis.
    F : numpy ndarray
        Array containing the values of the fluctuations in each window.
    nRho : numpy ndarray
        Array of window's sizes used for the computation of `rho`.
    rho : numpy ndarray
        Array containing the cross-correlation index in each window.
    nThr : numpy ndarray
        Array of window's sizes used for the computation of `rho` thresholds.
    confUp : numpy ndarray
        Array containing the first confidence interval in each window.
    confDown : numpy ndarray
        Array containing the second confidence interval in each window.
    isComputed : bool
        Boolean value to know if `F` has been computed in order to prevent the
        computation of other functions that need `F`.
    """

    def __init__(self, tsVec1=[], tsVec2=[]):
    	pass

    def computeFlucVec(self, winSizes, polOrd=1, absVals=True, overlap=False, revSeg=False):
        """Computation of the fluctuations in each window.

        Parameters
        ----------
        winSizes : numpy ndarray
            Array of window's sizes.
        polOrd : int, optional
            Order of the polynomial to be fitted in each window (default : 1).
        absVals : bool, optional
            If True, the computation of `F` is performed using the abolute values of
            the fluctuations of both `tsVec1` and `tsVec2` (default : True).
        overlap : bool, optional
            If True, computes `F` using overlapping segments (default : False).
        revSeg : bool, optional
            If True, the computation of `F` is repeated starting from the end of
            the time series, ignored if `overlap` is True (default : False).

        Returns
        -------
        numpy ndarray
            Array `n` of window's sizes.
        numpy ndarray
            Array `F` containing the values of the fluctuations in each window.
        """
        return 0

    def fitFlucVec(self, nStart=-999, nEnd=-999, logBase=np.e, verbose=False):
        """Fit of the fluctuations values.

        Parameters
        ----------
        nStart : int, optional
            Size of the smaller window used to fit `F` (default : first value of `n`).
        nEnd : int, optional
            Size of the bigger window used to fit `F` (default : last value of `n`).
        logBase : float, optional
            Base of the logarithm for the log-log fit of `n` vs `F` (default : e).
        verbose : bool, optional
            Verbosity (default : False).

        Returns
        -------
        float
            Slope of the fit.
        float
            Intercept of the fit.
        """
        return 0

    def multiFitFlucVec(self, limitsList, logBase=np.e, verbose=False):
        """Fit of the fluctuations values in different intervals at the same time.

        Parameters
        ----------
        limitsList : numpy ndarray
            kx2 array with the sizes of k starting and ending windows used to fit `F`.
        logBase : float, optional
            Base of the logarithm for the log-log fit of `n` vs `F` (default : e).
        verbose : bool, optional
            Verbosity (default : False).

        Returns
        -------
        numpy ndarray
            Slopes of the fits.
        numpy ndarray
            Intercepts of the fits.
        """
        return 0

    def computeRho(self, winSizes, polOrd=1, verbose=False, overlap=False, revSeg=False):
        """Computation of the cross-correlation index in each window.

        Parameters
        ----------
        winSizes : numpy ndarray
            Array of window's sizes.
        polOrd : int, optional
            Order of the polynomial to be fitted in each window (default : 1).
        verbose : bool, optional
            Verbosity (default : False).
        overlap : bool, optional
            If True, computes `F` using overlapping segments (default : False).
        revSeg : bool, optional
            If True, the computation of `F` is repeated starting from the end of
            the time series, ignored if `overlap` is True (default : False).

        Returns
        -------
        numpy ndarray
            Array of window's sizes.
        numpy ndarray
            Array containing the cross-correlation index.
        """
        return 0

    def rhoThresholds(self, winSizes, nSim, confLvl, polOrd=1, verbose=False):
        """Computation of the cross-correlation index's confidence levels in each window.

        Parameters
        ----------
        L : int
            Size of the random time series used to evaluate confidence levels.
        winSizes : numpy ndarray
            Array of window's sizes.
        nSim : int
            Number of times the cross-correlation index between two random time series
            is computed in order to evaluate the confidence levels.
        confLvl : float
            Confidence level.
        polOrd : int, optional
            Order of the polynomial to be fitted in each window (default : 1).
        verbose : bool, optional
            Verbosity (default : False).

        Returns
        -------
        numpy ndarray
            Array of window's sizes.
        numpy ndarray
            Array containing the first confidence interval.
        numpy ndarray
            Array containing the second confidence interval.
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

