#    ht.pyx - ht algorithm of fathon package
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

class HT:
    """Time-dependent local Hurst exponent class.
    
    Parameters
    ----------
    tsVec : iterable
        Time series used for the analysis.
    """

    def __init__(self, tsVec):
        pass

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
        return 0
