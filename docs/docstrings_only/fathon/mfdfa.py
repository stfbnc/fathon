#    mfdfa.pyx - mfdfa algorithm of fathon package
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
class MFDFA:
    """MultiFractal Detrended Fluctuation Analysis class.

    Parameters
    ----------
    n : numpy ndarray
        Array of window's sizes used for the computation.
    tsVec : iterable
        Time series used for the analysis.
    F : numpy ndarray
        Array containing the values of the fluctuations in every window.
    list_H : numpy ndarray
        Array containing the values of the slope of the fit at every q-order.
    q_list : numpy ndarray
        Array containing the values of the q-orders.
    nStep : int
        Value of the step between two consecutive window's sizes in `n`.
    isComputed : bool
        Boolean value to know if `F` has been computed in order to prevent the computation of other functions that need `F`.
    """

    def __init__(self, tsVec):
    	pass

    def computeFlucVec(self, nMin, q_list, nMax=-999, polOrd=1, nStep=1, revSeg=False):
        """Computation of the fluctuations in every window for every q-order.

        Parameters
        ----------
        nMin : int
            Size of the smaller window used to compute `F`.
        q_list : float or iterable or numpy ndarray
            List of q-orders used to compute `F`.
        nMax : int, optional
            Size of the bigger window used to compute `F` (default : len(`tsVec`)/4)).
        polOrd : int, optional
            Order of the polynomial to be fitted in every window (default : 1).
        nStep : int, optional
            Step between two consecutive window's sizes (default : 1).
        revSeg : bool, optional
            If True, the computation of `F` is repeated starting from the end of the time series (default : False).

        Returns
        -------
        numpy ndarray
            Array `n` of window's sizes.
        numpy ndarray
            qxn array `F` containing the values of the fluctuations in every window for every q-order.
        """
        return 0

    def fitFlucVec(self, n_start=-999, n_end=-999):
        """Fit of the fluctuations values.

        Parameters
        ----------
        n_start : int, optional
            Size of the smaller window used to fit `F` at every q-order (default : first value of `n`).
        n_end : int, optional
            Size of the bigger window used to fit `F` at every q-order (default : last value of `n`).

        Returns
        -------
        numpy ndarray
            Slope of the fit for every q-order.
        numpy ndarray
            Intercept of the fit for every q-order.
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

