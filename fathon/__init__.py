#    __init__.py - init for fathon package
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

import platform
import os
init_dir = os.path.dirname(os.path.abspath(__file__))
if platform.system() == "Linux":
    cmd = "cp " + init_dir + "/3rd_party/gsl/lib/*.so* " + init_dir + "/../fathon.libs"
    os.system(cmd)
    
from .dfa import DFA
from .mfdfa import MFDFA
from .dcca import DCCA
from .ht import HT

"""A Python package for detrended fluctuation analysis (DFA)
    and related algorithms.
"""

__version__ = '1.0'
__author__ = 'Stefano Bianchi'
__git_repo__ = 'https://github.com/stfbnc/fathon'
