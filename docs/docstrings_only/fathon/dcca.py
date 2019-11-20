#    dcca.pyx - dcca algorithm of fathon package
#    Copyright (C) 2019  Stefano Bianchi
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
        Array containing the values of the fluctuations in every window.
    nStep : int
        Value of the step between two consecutive window's sizes in `n`.
    isComputed : bool
        Boolean value to know if `F` has been computed in order to prevent the computation of other functions that need `F`.
    """

    def __init__(self, tsVec1=[], tsVec2=[]):
        pass

    def computeFlucVec(self, nMin, nMax=-999, polOrd=1, nStep=1, absVals=True):
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
        return 0

    def fitFlucVec(self, n_start=-999, n_end=-999):
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
        return 0

    def multiFitFlucVec(self, limits_list):
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
        return 0

    def computeRho(self, nMin, nMax=-999, polOrd=1, nStep=1):
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
        return 0

    def rhoThresholds(self, L, nMin, nMax, nSim, confLvl, polOrd=1, nStep=1):
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
        return 0
