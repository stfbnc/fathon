#    dcca.pyx - dcca algorithm of fathon package
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
import warnings

cdef extern from "cLoops.h" nogil:
    double flucDCCAAbsCompute(double *y1, double *y2, int curr_win_size, int N, int pol_ord)
    double flucDCCANoAbsCompute(double *y1, double *y2, int curr_win_size, int N, int pol_ord)

cdef class DCCA:
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
        Array containing the values of the fluctuations in every window.
    nStep : int
        Value of the step between two consecutive window's sizes in `n`.
    isComputed : bool
        Boolean value to know if `F` has been computed in order to prevent the computation of other functions that need `F`.
    """

    cdef:
        np.ndarray n
        np.ndarray tsVec1, tsVec2, F
        int nStep
        bint isComputed

    def __init__(self, tsVec1=[], tsVec2=[]):
        if len(tsVec1) != 0 and len(tsVec2) != 0:
            self.tsVec1 = np.array(tsVec1, dtype=float)
            self.tsVec1 = self.tsVec1[~np.isnan(self.tsVec1)]
            self.tsVec2 = np.array(tsVec2, dtype=float)
            self.tsVec2 = self.tsVec2[~np.isnan(self.tsVec2)]
            if len(self.tsVec1) != len(self.tsVec2):
                warnings.warn("Warning: Input vectors have different length. The longest vector has been reduced to the size of the shortest one.")
                self.tsVec1 = self.tsVec1[0:np.min([len(self.tsVec1), len(self.tsVec2)])]
                self.tsVec2 = self.tsVec2[0:np.min([len(self.tsVec1), len(self.tsVec2)])]
        self.isComputed = False

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cdef cy_flucCompute(self, np.ndarray[np.float64_t, ndim=1, mode='c'] vects1, np.ndarray[np.float64_t, ndim=1, mode='c'] vects2,
                                np.ndarray[int, ndim=1, mode='c'] vecn, np.ndarray[np.float64_t, ndim=1, mode='c'] vecf, int polOrd, bint absVals):
        cdef int nLen, tsLen
        cdef Py_ssize_t i

        nLen = len(vecn)
        tsLen = len(vects1)
        with nogil:
            if absVals:
                for i in range(nLen):
                    vecf[i] = flucDCCAAbsCompute(&vects1[0], &vects2[0], vecn[i], tsLen, polOrd)
            else:
                for i in range(nLen):
                    vecf[i] = flucDCCANoAbsCompute(&vects1[0], &vects2[0], vecn[i], tsLen, polOrd)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cpdef computeFlucVec(self, int nMin, int nMax=-999, int polOrd=1, int nStep=1, bint absVals=True):
        """Computation of the fluctuations in every window.

        Parameters
        ----------
        nMin : int
            Size of the smaller window used to compute `F`.
        nMax : int, optional
            Size of the bigger window used to compute `F` (default : len(`tsVec`)/4)).
        polOrd : int, optional
            Order of the polynomial to be fitted in every window (default : 1).
        nStep : int, optional
            Step between two consecutive window's sizes (default : 1).
        absVals : bool, optional
            If True, the computation of `F` is performed using the abolute values of the fluctuations of both `tsVec1` and `tsVec2` (default : True).

        Returns
        -------
        numpy ndarray
            Array `n` of window's sizes.
        numpy ndarray
            Array `F` containing the values of the fluctuations in every window.
        """
        cdef int tsLen = len(self.tsVec1)

        if polOrd < 1:
            raise SystemExit('Error: Polynomial order must be greater than 0.')
        if nStep < 1:
            raise SystemExit('Error: Step for scales must be greater than 0.')
        if nMax == -999:
            nMax = int(tsLen/4)
        if nMax < 3 or nMin < 3:
            raise SystemExit('Error: Variable nMin and nMax must be at least equal to 3.')
        if nMax <= nMin:
            raise SystemExit('Error: Variable nMax must be greater than variable nMin.')
        if nMax > tsLen:
            raise SystemExit('Error: Variable nMax must be less than the input vector length.')
        if nMin < (polOrd+2):
            raise SystemExit('Error: Variable nMin must be at least equal to {}.'.format(polOrd+2))

        self.nStep = nStep
        self.n = np.arange(nMin, nMax+1, nStep, dtype=ctypes.c_int)
        self.F = np.zeros((len(self.n), ), dtype=ctypes.c_double)
        self.cy_flucCompute(np.array(self.tsVec1, dtype=ctypes.c_double), np.array(self.tsVec2, dtype=ctypes.c_double), self.n, self.F, polOrd, absVals)
        self.isComputed = True
        return self.n, self.F

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cdef computeFlucVecSameTs(self, np.ndarray[np.float64_t, ndim=1, mode='c'] vec, int nMin, int nMax, int polOrd):
        cdef int nLen, tsLen = len(vec)
        cdef Py_ssize_t i
        cdef np.ndarray[np.float64_t, ndim=1, mode='c'] F_same
        cdef np.ndarray[int, ndim=1, mode='c'] vecn

        vecn = np.arange(nMin, nMax+1, self.nStep, dtype=ctypes.c_int)
        nLen = len(vecn)
        F_same = np.zeros((nLen, ), dtype=ctypes.c_double)
        with nogil:
            for i in range(nLen):
                F_same[i] = flucDCCAAbsCompute(&vec[0], &vec[0], vecn[i], tsLen, polOrd)
        return F_same

    @cython.boundscheck(False)
    @cython.nonecheck(False)
    cpdef fitFlucVec(self, int n_start=-999, int n_end=-999):
        """Fit of the fluctuations values.

        Parameters
        ----------
        n_start : int, optional
            Size of the smaller window used to fit `F` (default : first value of `n`).
        n_end : int, optional
            Size of the bigger window used to fit `F` (default : last value of `n`).

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
            if n_start == -999:
                n_start = self.n[0]
            if n_end == -999:
                n_end = self.n[-1]
            if n_start > n_end:
                raise SystemExit('Error: Variable n_end must be greater than variable n_start.')
            if (n_start < self.n[0]) or (n_end > self.n[-1]):
                raise SystemExit('Error: Fit limits must be included in interval [{}, {}].'.format(self.n[0], self.n[-1]))
            if (n_start not in self.n) or (n_end not in self.n):
                raise SystemExit('Error: Fit limits must be included in the n vector.')

            start = int((n_start - self.n[0]) / self.nStep)
            end = int((n_end - self.n[0]) / self.nStep)
            log_fit = np.polyfit(np.log(self.n[start:end+1]) , np.log(self.F[start:end+1]), 1)
            print('Fit limits: [{}, {}]'.format(n_start, n_end))
            print('Fit result: H_intercept = {:.2f}, H = {:.2f}'.format(log_fit[1], log_fit[0]))
            return log_fit[0], log_fit[1]
        else:
            raise SystemExit('Error: Fluctuations vector has not been computed yet.')

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cpdef multiFitFlucVec(self, np.ndarray[np.int_t, ndim=2, mode='c'] limits_list):
        """Fit of the fluctuations values in different intervals at the same time.

        Parameters
        ----------
        limits_list : numpy ndarray
            kx2 array with the sizes of k starting and ending windows used to fit `F`.

        Returns
        -------
        numpy ndarray
            Slopes of the fits.
        numpy ndarray
            Intercepts of the fits.
        """
        cdef Py_ssize_t i
        cdef int limLen = len(limits_list)
        cdef np.ndarray[np.float64_t, ndim=1, mode='c'] list_H_intercept, list_H

        if self.isComputed:
            list_H = np.zeros((limLen, ), dtype=float)
            list_H_intercept = np.zeros((limLen, ), dtype=float)
            for i in range(limLen):
                print('----------')
                list_H[i], list_H_intercept[i] = self.fitFlucVec(n_start=limits_list[i][0], n_end=limits_list[i][1])
            print('----------')
            return list_H, list_H_intercept
        else:
            raise SystemExit('Error: Fluctuations vector has not been computed yet.')

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cpdef computeRho(self, int nMin, int nMax=-999, int polOrd=1, int nStep=1):
        """Computation of the cross-correlation index in every window.

        Parameters
        ----------
        nMin : int
            Size of the smaller window used to compute the cross-correlation index.
        nMax : int, optional
            Size of the bigger window used to compute the cross-correlation index (default : len(`tsVec`)/4)).
        polOrd : int, optional
            Order of the polynomial to be fitted in every window (default : 1).
        nStep : int, optional
            Step between two consecutive window's sizes (default : 1).

        Returns
        -------
        numpy ndarray
            Array of window's sizes.
        numpy ndarray
            Array containing the cross-correlation index.
        """
        cdef Py_ssize_t i
        cdef int nLen, tsLen = len(self.tsVec1)
        cdef np.ndarray[int, ndim=1, mode='c'] nn
        cdef np.ndarray[np.float64_t, ndim=1, mode='c'] Fxy, Fxx, Fyy, rho

        if polOrd < 1:
            raise SystemExit('Error: Polynomial order must be greater than 0.')
        if nStep < 1:
            raise SystemExit('Error: Step for scales must be greater than 0.')
        if nMax == -999:
            nMax = int(tsLen/4)
        if nMax < 3 or nMin < 3:
            raise SystemExit('Error: Variable nMin and nMax must be at least equal to 3.')
        if nMax <= nMin:
            raise SystemExit('Error: Variable nMax must be greater than variable nMin.')
        if nMax > tsLen:
            raise SystemExit('Error: Variable nMax must be less than the input vector length.')
        if nMin < (polOrd+2):
            raise SystemExit('Error: Variable nMin must be at least equal to {}.'.format(polOrd+2))

        self.nStep = nStep
        nn, Fxy = self.computeFlucVec(nMin, nMax=nMax, polOrd=polOrd, nStep=self.nStep, absVals=False)
        print('DCCA between series 1 and 2 computed.')
        Fxx = self.computeFlucVecSameTs(self.tsVec1, nMin, nMax=nMax, polOrd=polOrd)
        print('DCCA between series 1 and 1 computed.')
        Fyy = self.computeFlucVecSameTs(self.tsVec2, nMin, nMax=nMax, polOrd=polOrd)
        print('DCCA between series 2 and 2 computed.')

        nLen = len(nn)
        rho = np.zeros((nLen, ), dtype=float)
        with nogil:
            for i in range(nLen):
                rho[i] = Fxy[i] / (Fxx[i] * Fyy[i])
        return nn, rho

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cpdef rhoThresholds(self, int L, int nMin, int nMax, int nSim, double confLvl, int polOrd=1, int nStep=1):
        """Computation of the cross-correlation index's confidence levels in every window.

        Parameters
        ----------
        L : int
            Size of the random time series used to evaluate confidence levels.
        nMin : int
            Size of the smaller window used to compute the cross-correlation index.
        nMax : int
            Size of the bigger window used to compute the cross-correlation index.
        nSim : int
            Number of times the cross-correlation index between two random time series is computed in order to evaluate the confidence levels.
        confLvl : float
            Confidence level.
        polOrd : int, optional
            Order of the polynomial to be fitted in every window (default : 1).
        nStep : int, optional
            Step between two consecutive window's sizes (default : 1).

        Returns
        -------
        numpy ndarray
            Array of window's sizes.
        numpy ndarray
            Array containing the first confidence interval.
        numpy ndarray
            Array containing the second confidence interval.
        """
        cdef np.ndarray[np.float64_t, ndim=2, mode='c'] rho_all
        cdef np.ndarray[np.float64_t, ndim=1, mode='c'] ran1, ran2, vecfx, vecfy, vecfxy, up_lim, down_lim
        cdef np.ndarray[int, ndim=1, mode='c'] vecn
        cdef int nLen

        if polOrd < 1:
            raise SystemExit('Error: Polynomial order must be greater than 0.')
        if nStep < 1:
            raise SystemExit('Error: Step for scales must be greater than 0.')
        if nMax < 3 or nMin < 3:
            raise SystemExit('Error: Variable nMin and nMax must be at least equal to 3.')
        if nMax <= nMin:
            raise SystemExit('Error: Variable nMax must be greater than variable nMin.')
        if nMin < (polOrd+2):
            raise SystemExit('Error: Variable nMin must be at least equal to {}.'.format(polOrd+2))
        if nSim < 1:
            raise SystemExit('Error: Number of simulations must be greater than 0.')
        if confLvl < 0 or confLvl > 1:
            raise SystemExit('Error: Confidence level must be incliuded in the interval [0,1].')

        vecn = np.arange(nMin, nMax+1, nStep, dtype=ctypes.c_int)
        nLen = len(vecn)
        up_lim = np.zeros((nLen, ), dtype=ctypes.c_double)
        down_lim = np.zeros((nLen, ), dtype=ctypes.c_double)
        rho_all = np.zeros((nSim, nLen), dtype=ctypes.c_double)
        for i in range(nSim):
            print('Simulation number {}'.format(i+1))
            ran1 = np.random.randn(L)
            ran1 = np.cumsum(ran1-np.mean(ran1))
            ran2 = np.random.randn(L)
            ran2 = np.cumsum(ran2-np.mean(ran2))
            vecfx = np.zeros((nLen, ), dtype=ctypes.c_double)
            vecfy = np.zeros((nLen, ), dtype=ctypes.c_double)
            vecfxy = np.zeros((nLen, ), dtype=ctypes.c_double)
            self.cy_flucCompute(ran1, ran2, vecn, vecfxy, polOrd, False)
            self.cy_flucCompute(ran1, ran1, vecn, vecfx, polOrd, False)
            self.cy_flucCompute(ran2, ran2, vecn, vecfy, polOrd, False)
            for j in range(nLen):
                rho_all[i, j] = vecfxy[j] / (vecfx[j] * vecfy[j])
        up_lim = np.quantile(rho_all, confLvl, axis=0)
        down_lim = np.quantile(rho_all, 1-confLvl, axis=0)
        return vecn, up_lim, down_lim
