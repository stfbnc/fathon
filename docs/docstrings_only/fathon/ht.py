#    ht.pyx - ht algorithm of fathon package
#    Copyright (C) 2019-2020  Stefano Bianchi
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
from . import mfdfa
from . import fathonUtils as fu
class HT:
    """Time-dependent local Hurst exponent class.
    
    Parameters
    ----------
    tsVec : iterable
        Time series used for the analysis.
    ht : numpy ndarray
        Time-dependent local Hurst exponent.
    """

    def __init__(self, tsVec):
    	pass

    def computeHt(self, scales, polOrd=1, mfdfaPolOrd=1, q0Fit=[], verbose=False):
        """Computation of the time-dependent local Hurst exponent at every scale, using Ihlen's approach.
        
        Parameters
        ----------
        scales : int or iterable or numpy ndarray
            Window's sizes used for the computation of the time-dependent local Hurst exponent.
        polOrd : int, optional
            Order of the polynomial to be fitted in every window (default : 1).
        mfdfaPolOrd : int, optional
            Order of the polynomial to be fitted to MFDFA's fluctuations at q = 0 (default : 1).
        q0Fit : iterable or numpy ndarray of floats, optional
            MFDFA's Hurst exponent at order q = 0 and the corresponding intercept of the fit, [hq0, hq0_intercept]. These values must come from a log-log fit, with the log base equal to e. If not empty, it will be directly used to compute the time-dependent local Hurst exponent, ignoring `mfdfaPolOrd` value (default : []).
        verbose : bool, optional
            Verbosity (default : False).

        Returns
        -------
        numpy ndarray
            Time-dependent local Hurst exponent.
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

