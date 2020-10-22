#    fathonUtils.py - utils functions for fathon package
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
import pickle

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

def linRangeByStep(start, end, step=1):
    """Array of linearly separated elements.
    
    Parameters
    ----------
    start : int
        Smallest element.
    end : int
        Biggest element.
    step : int
        Step between two consecutive elements (default : 1).

    Returns
    -------
    numpy ndarray
        Array of linearly separated elements.
    """
    return np.arange(start, end + 1, step, dtype=int)

def linRangeByCount(start, end, count=-1):
    """Array of linearly separated elements.
    
    Parameters
    ----------
    start : int
        Smallest element.
    end : int
        Biggest element.
    count : int
        Number of elements (default : `end` - `start` + 1).

    Returns
    -------
    numpy ndarray
        Array of linearly separated elements.
    """
    if count == -1 or count > (end - start + 1):
        count = end - start + 1
        
    return np.linspace(start, end, count, endpoint=True, dtype=int)

def powRangeByStep(start, end, step=1, base=2):
    """Array of elements given by `base` raised to linearly separated exponents.
    
    Parameters
    ----------
    start : int
        Smallest element.
    end : int
        Biggest element.
    step : int
        Step between two consecutive exponents (default : 1).
    base : int
        Base of the exponential (default : 2).

    Returns
    -------
    numpy ndarray
        Array of elements given by `base` raised to linearly separated exponents.
    """
    exponents = linRangeByStep(start, end, step=step)
    
    return np.power(base, exponents, dtype=int)

def powRangeByCount(start, end, count=-1, base=2):
    """Array of elements given by `base` raised to linearly separated exponents.
    
    Parameters
    ----------
    start : int
        Smallest element.
    end : int
        Biggest element.
    count : int
        Number of elements (default : `end` - `start` + 1).
    base : int
        Base of the exponential (default : 2).

    Returns
    -------
    numpy ndarray
        Array of elements given by `base` raised to linearly separated exponents.
    """
    exponents = linRangeByCount(start, end, count=count)
        
    return np.power(base, exponents, dtype=int)

def getObjectMember(fileName, memberName):
    """Return member of a previously saved object. Member's name is the same of the object's member it refers to. Member `isComputed` has no practical use and cannot be retrieved.

    Parameters
    ----------
    fileName : str
        Path to a previously saved `.fathon` file.
    memberName : str
        Desired member's name.

    Returns
    -------
    numpy ndarray
        Object's member.
    """
    if fileName.split('.')[-1] != 'fathon':
        raise ValueError('Error: Not recognized extension.')

    f = open(fileName, 'rb')
    data = pickle.load(f)
    f.close()
    if memberName not in data.keys() or memberName == 'isComputed':
        raise ValueError('Error: Member not present.')
    else:
        ret = data[memberName]
        if len(ret) != 0:
            ret = np.array(ret, dtype=type(ret[0]))

    return ret

