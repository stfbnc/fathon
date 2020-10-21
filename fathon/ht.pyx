#    ht.pyx - ht algorithm of fathon package
#    Copyright (C) 2019-2020  Stefano Bianchi
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

#cython: language_level=3

import numpy as np
cimport numpy as np
cimport cython
from cython.parallel import prange
import ctypes
import pickle
from . import mfdfa
from . import fathonUtils as fu
	
cdef extern from "cLoops.h" nogil:
    double HTCompute(double *y, int scale, int N, int pol_ord, int v)

cdef class HT:
    """Time-dependent local Hurst exponent class.
    
    Parameters
    ----------
    tsVec : iterable
        Time series used for the analysis.
    ht : numpy ndarray
        Time-dependent local Hurst exponent.
    """

    cdef:
        np.ndarray tsVec, ht

    def __init__(self, tsVec):
        if isinstance(tsVec, str):
            if len(tsVec.split('.')) > 1 and tsVec.split('.')[-1] == 'fathon':
                f = open(tsVec, 'rb')
                data = pickle.load(f)
                f.close()
                if data['kind'] != 'ht':
                    raise ValueError('Error: Loaded object is not a HT object.')
                else:
                    self.tsVec = np.array(data['tsVec'], dtype=float)
                    self.ht = np.array(data['ht'], dtype=ctypes.c_double)
            else:
                raise ValueError('Error: Not recognized extension.')
        else:
            self.tsVec = np.array(tsVec, dtype=float)
            self.tsVec = self.tsVec[~np.isnan(self.tsVec)]
		
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cdef cy_computeHt(self, np.ndarray[int, ndim=1, mode='c'] scales, int polOrd, int mfdfaPolOrd, np.ndarray[np.float64_t, ndim=1, mode='c'] q0Fit, bint verbose):
        cdef int htRowLen, tsLen, scale
        cdef Py_ssize_t i, j
        cdef double H0, H0_intercept
        cdef np.ndarray[np.float64_t, ndim=1, mode='c'] vects, vecht
        
        tsLen = len(self.tsVec)
        htRowLen = tsLen - min(scales) + 1
        vects = np.array(self.tsVec, dtype=ctypes.c_double)
        vecht = np.zeros((htRowLen * len(scales), ), dtype=ctypes.c_double)
        
        if len(q0Fit) == 0:
            pymfdfa = mfdfa.MFDFA(self.tsVec)
            _, _ = pymfdfa.computeFlucVec(fu.linRangeByCount(10, int(tsLen / 4), count=20), 0.0, revSeg=True, polOrd=mfdfaPolOrd)
            H0, H0_intercept = pymfdfa.fitFlucVec(verbose=verbose)
        else:
            if verbose:
                print('Variable q0Fit assigned, variable mfdfaPolOrd will be ignored.')
            H0 = q0Fit[0]
            H0_intercept = q0Fit[1]
        
        if verbose:
            print('-----')
        for i, scale in enumerate(scales):
            if verbose:
                print('scale = {}'.format(scale))
                
            for j in prange(tsLen - scale + 1, nogil=True):
               vecht[i*htRowLen+j] = HTCompute(&vects[0], scale, tsLen, polOrd, j)
               
            if verbose:
                print('-----')
                
            for j in range(htRowLen):
                if vecht[i*htRowLen+j] != 0.0:
                    vecht[i*htRowLen+j] = (H0_intercept + H0 * np.log(scale) - np.log(vecht[i*htRowLen+j])) / (np.log(tsLen - scale + 1) - np.log(scale)) + H0
                    
        return np.reshape(vecht, (len(scales), htRowLen))
		
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
        if polOrd < 1:
            raise ValueError('Error: Polynomial order must be greater than 0.')
        if mfdfaPolOrd < 1:
            raise ValueError('Error: Polynomial order must be greater than 0.')
        if len(q0Fit) == 1 or len(q0Fit) > 2:
            raise ValueError('Error: List q0Fit must have length 2.')
        if len(q0Fit) == 2 and q0Fit[0] <= 0.0:
            raise ValueError('Error: First element of q0Fit must be greater than 0.')

        if isinstance(scales, int):
            if scales < 3:
                raise ValueError('Error: Every scale must be at least equal to 3.')
            else:
                scales = np.array([scales], dtype=ctypes.c_int)
        elif isinstance(scales, list) or isinstance(scales, np.ndarray):
            for scale in scales:
                if scale < 3:
                    raise ValueError('Error: Every scale must be at least equal to 3.')
            scales = np.array(scales, dtype=ctypes.c_int)
        else:
            raise ValueError('Error: scales type is {}. Expected int, list, or numpy array.'.format(type(scales)))
         
        q0Fit = np.array(q0Fit, dtype=ctypes.c_double)
        self.ht = self.cy_computeHt(scales, polOrd, mfdfaPolOrd, q0Fit, verbose)
        
        return self.ht

    def saveObject(self, outFileName):
        """Save current object state to binary file.
        
        Parameters
        ----------
        outFileName : str
            Output binary file. `.fathon` extension will be appended to the file name.
        """
        saveDict = {}
        saveDict['kind'] = 'ht'
        saveDict['tsVec'] = self.tsVec.tolist()
        try:
            saveDict['ht'] = self.ht.tolist()
        except:
            saveDict['ht'] = []

        f = open(outFileName + '.fathon', 'wb')
        pickle.dump(saveDict, f)
        f.close()
