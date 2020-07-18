#    tsHelper.py - helper for fathon package
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

import numpy as np
import ctypes

def subtractMean(vec):
	"""Subtracts mean of a vector.
	
	Parameters
	----------
	vec : iterable
		Python iterable whose mean will be subtracted.
	
	Returns
	-------
	numpy ndarray
		`vec` with its mean subtracted.
	"""
	return np.array(vec - np.nanmean(vec), dtype=float)
	
def toAggregated(vec):
	"""Subtracts mean of a vector and computes the cumulative sum.
	
	Parameters
	----------
	vec : iterable
		The vector to be integrated.
	
	Returns
	-------
	numpy ndarray
		Cumulative sum of `vec` after subtraction of its mean.
	"""
	return np.array(np.nancumsum(vec - np.nanmean(vec)), dtype=float)

def windowsVec(start, end, step):
    """Linear range of window's sizes
    
    Parameters
    ----------
    start : int
        Smallest size
    end : int
        Biggest size
    """
        
# o scale separate da step o numero di scale predefinito