#    dfa.pyx - dfa algorithm of fathon package
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
import ctypes
import pickle

cdef extern from "cLoops.h" nogil:
    double flucDFAForwCompute(double *y, int curr_win_size, int N, int pol_ord)
    double flucDFAForwBackwCompute(double *y, int curr_win_size, int N, int pol_ord)

cdef class DFA:
    """Detrended Fluctuation Analysis class.

    Parameters
    ----------
    n : numpy ndarray
        Array of window's sizes used for the computation.
    tsVec : iterable
        Time series used for the analysis.
    F : numpy ndarray
        Array containing the values of the fluctuations in every window.
    isComputed : bool
        Boolean value to know if `F` has been computed in order to prevent the computation of other functions that need `F`.
    """

    cdef:
        np.ndarray n
        np.ndarray tsVec, F
        bint isComputed

    def __init__(self, tsVec):
        if isinstance(tsVec, str):
            if len(tsVec.split('.')) > 1 and tsVec.split('.')[-1] == 'fathon':
                f = open(tsVec, 'rb')
                data = pickle.load(f)
                f.close()
                if data['kind'] != 'dfa':
                    raise ValueError('Error: Loaded object is not a DFA object.')
                else:
                    self.tsVec = np.array(data['tsVec'], dtype=float)
                    self.n = np.array(data['n'], dtype=ctypes.c_int)
                    self.F = np.array(data['F'], dtype=ctypes.c_double)
                    self.isComputed = data['isComputed']
            else:
                raise ValueError('Error: Not recognized extension.')
        else:
            self.tsVec = np.array(tsVec, dtype=float)
            self.tsVec = self.tsVec[~np.isnan(self.tsVec)]
            self.isComputed = False

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cdef cy_flucCompute(self, np.ndarray[np.float64_t, ndim=1, mode='c'] vects, np.ndarray[int, ndim=1, mode='c'] vecn,
                        np.ndarray[np.float64_t, ndim=1, mode='c'] vecf, int polOrd, bint revSeg):
        cdef int nLen, tsLen
        cdef Py_ssize_t i

        nLen = len(vecn)
        tsLen = len(vects)
        with nogil:
            if revSeg:
                for i in range(nLen):
                    vecf[i] = flucDFAForwBackwCompute(&vects[0], vecn[i], tsLen, polOrd)
            else:
                for i in range(nLen):
                    vecf[i] = flucDFAForwCompute(&vects[0], vecn[i], tsLen, polOrd)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cpdef computeFlucVec(self, np.ndarray[np.int64_t, ndim=1, mode='c'] winSizes, int polOrd=1, bint revSeg=False):
        """Computation of the fluctuations in every window.

        Parameters
        ----------
        winSizes : numpy ndarray
            Array of window's sizes.
        polOrd : int, optional
            Order of the polynomial to be fitted in every window (default : 1).
        revSeg : bool, optional
            If True, the computation of `F` is repeated starting from the end of the time series (default : False).

        Returns
        -------
        numpy ndarray
            Array `n` of window's sizes.
        numpy ndarray
            Array `F` containing the values of the fluctuations in every window.
        """
        cdef int tsLen = len(self.tsVec)

        if polOrd < 1:
            raise ValueError('Error: Polynomial order must be greater than 0.')
        if winSizes[len(winSizes)-1] <= winSizes[0]:
            raise ValueError('Error: `winSizes[-1]` must be greater than variable `winSizes[0]`.')
        if winSizes[len(winSizes)-1] > tsLen:
            raise ValueError('Error: `winSizes[-1]` must be smaller than the input vector length.')
        if winSizes[0] < (polOrd + 2):
            raise ValueError('Error: `winSizes[0]` must be at least equal to {}.'.format(polOrd + 2))

        self.n = np.array(winSizes, dtype=ctypes.c_int)
        self.F = np.zeros((len(self.n), ), dtype=ctypes.c_double)
        self.cy_flucCompute(np.array(self.tsVec, dtype=ctypes.c_double), self.n, self.F, polOrd, revSeg)
        self.isComputed = True
        
        return self.n, self.F

    @cython.boundscheck(False)
    @cython.nonecheck(False)
    cpdef fitFlucVec(self, int nStart=-999, int nEnd=-999, float logBase=np.e, bint verbose=False):
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
        cdef int start, end
        cdef np.ndarray[np.float64_t, ndim=1, mode='c'] log_fit

        if self.isComputed:
            if nStart == -999:
                nStart = self.n[0]
            if nEnd == -999:
                nEnd = self.n[-1]
            if nStart > nEnd:
                raise ValueError('Error: Variable nEnd must be greater than variable nStart.')
            if (nStart < self.n[0]) or (nEnd > self.n[-1]):
                raise ValueError('Error: Fit limits must be included in interval [{}, {}].'.format(self.n[0], self.n[-1]))
            if (nStart not in self.n) or (nEnd not in self.n):
                raise ValueError('Error: Fit limits must be included in the window\'s sizes vector.')

            start = np.where(self.n==nStart)[0][0]
            end = np.where(self.n==nEnd)[0][0]
            log_fit = np.polyfit(np.log(self.n[start:end+1]) / np.log(logBase), np.log(self.F[start:end+1]) / np.log(logBase), 1)
            
            if verbose:
                print('Fit limits: [{}, {}]'.format(nStart, nEnd))
                print('Fit result: H intercept = {:.2f}, H = {:.2f}'.format(log_fit[1], log_fit[0]))
                
            return log_fit[0], log_fit[1]
        else:
            print('Nothing to fit, fluctuations vector has not been computed yet.')

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cpdef multiFitFlucVec(self, np.ndarray[np.int_t, ndim=2, mode='c'] limitsList, float logBase=np.e, bint verbose=False):
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
        cdef Py_ssize_t i
        cdef int limLen = len(limitsList)
        cdef np.ndarray[np.float64_t, ndim=1, mode='c'] list_H, list_H_intercept

        if self.isComputed:
            list_H = np.zeros((limLen, ), dtype=float)
            list_H_intercept = np.zeros((limLen, ), dtype=float)
            
            for i in range(limLen):
                if verbose:
                    print('----------')
                list_H[i], list_H_intercept[i] = self.fitFlucVec(nStart=limitsList[i][0], nEnd=limitsList[i][1], logBase=logBase, verbose=verbose)
                
            if verbose:
                print('----------')
                
            return list_H, list_H_intercept
        else:
            print('Nothing to fit, fluctuations vector has not been computed yet.')

    def saveObject(self, outFileName):
        """Save current object state to binary file.
        
        Parameters
        ----------
        outFileName : str
            Output binary file. `.fathon` extension will be appended to the file name.
        """
        saveDict = {}
        saveDict['kind'] = 'dfa'
        saveDict['tsVec'] = self.tsVec.tolist()
        try:
            saveDict['n'] = self.n.tolist()
        except:
            saveDict['n'] = []
        try:
            saveDict['F'] = self.F.tolist()
        except:
            saveDict['F'] = []
        saveDict['isComputed'] = self.isComputed

        f = open(outFileName + '.fathon', 'wb')
        pickle.dump(saveDict, f)
        f.close()

