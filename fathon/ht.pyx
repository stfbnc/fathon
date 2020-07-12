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
from . import mfdfa
	
cdef extern from "cLoops.h" nogil:
    double HTCompute(double *y, int scale, int N, int pol_ord, int v)

cdef class HT:
    """Time-dependent local Hurst exponent class.
    
    Parameters
    ----------
    tsVec : iterable
        Time series used for the analysis.
    """

    cdef:
        np.ndarray tsVec

    def __init__(self, tsVec):
        self.tsVec = np.array(tsVec, dtype=float)
        self.tsVec = self.tsVec[~np.isnan(self.tsVec)]
		
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cdef cy_computeHt(self, np.ndarray[int, ndim=1, mode='c'] scales, int polOrd, int mfdfaPolOrd):
        cdef int htRowLen, tsLen, scale, mfdfa_step
        cdef Py_ssize_t i, j
        cdef double H0, H0_intercept
        cdef np.ndarray[np.float64_t, ndim=1, mode='c'] vects, vecht
        
        tsLen = len(self.tsVec)
        htRowLen = tsLen - min(scales) + 1
        vects = np.array(self.tsVec, dtype=ctypes.c_double)
        vecht = np.zeros((htRowLen*len(scales), ), dtype=ctypes.c_double)
        mfdfa_step = int(tsLen/100) if tsLen > 100 else 1
        print('-----')
        for i, scale in enumerate(scales):
            for j in prange(tsLen-scale+1, nogil=True):
               vecht[i*htRowLen+j] = HTCompute(&vects[0], scale, tsLen, polOrd, j)
            pymfdfa = mfdfa.MFDFA(self.tsVec)
            print('scale = {}'.format(scale))
            _, _ = pymfdfa.computeFlucVec(10, 0.0, nMax=int(tsLen/4), nStep=mfdfa_step, revSeg=True, polOrd=mfdfaPolOrd)
            H0, H0_intercept = pymfdfa.fitFlucVec()
            print('-----')
            for j in range(htRowLen):
                if vecht[i*htRowLen+j] != 0.0:
                    vecht[i*htRowLen+j] = (H0_intercept+H0*np.log(scale)-np.log(vecht[i*htRowLen+j])) / (np.log(tsLen-scale+1)-np.log(scale)) + H0
        return np.reshape(vecht, (len(scales), htRowLen))
		
    def computeHt(self, scales, polOrd=1, mfdfaPolOrd=1):
        """Computation of the time-dependent local Hurst exponent at every scale.
        
        Parameters
        ----------
        scales : int or iterable or numpy ndarray
            Window's sizes used for the computation of the time-dependent local Hurst exponent.
        polOrd : int, optional
            Order of the polynomial to be fitted in every window (default : 1).
        mfdfaPolOrd : int, optional
            Order of the polynomial to be fitted to MFDFA's fluctuations at q = 0 (default : 1).

        Returns
        -------
        numpy ndarray
            Time-dependent local Hurst exponent.
        """
        if polOrd < 1:
            raise SystemExit('Error: Polynomial order must be greater than 0.')
        if mfdfaPolOrd < 1:
            raise SystemExit('Error: Polynomial order must be greater than 0.')

        if isinstance(scales, int):
            if scales < 3:
                raise SystemExit('Error: Every scale at least equal to 3.')
            else:
                scales = np.array([scales], dtype=ctypes.c_int)
        elif isinstance(scales, list) or isinstance(scales, np.ndarray):
            for scale in scales:
                if scale < 3:
                    raise SystemExit('Error: Every scale at least equal to 3.')
            scales = np.array(scales, dtype=ctypes.c_int)
        else:
            raise SystemExit('Error: scales type is {}. Expected int, list, or numpy array.'.format(type(scales)))
        ht = self.cy_computeHt(scales, polOrd, mfdfaPolOrd)
        return ht
