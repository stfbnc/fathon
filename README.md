# fathon

[![Build Status](https://travis-ci.org/stfbnc/fathon.svg?branch=master)](https://travis-ci.org/stfbnc/fathon) [![Documentation Status](https://readthedocs.org/projects/fathon/badge/?version=latest)](https://fathon.readthedocs.io/en/latest/?badge=latest) [![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-orange.svg)](https://opensource.org/licenses/) [![Issues](https://img.shields.io/github/issues-raw/stfbnc/fathon.svg?maxAge=25000)](https://github.com/stfbnc/fathon/issues) [![PyPI version](https://badge.fury.io/py/fathon.svg)](https://badge.fury.io/py/fathon)  [![GitHub stars](https://img.shields.io/github/stars/stfbnc/fathon.svg?style=social&label=Stars&style=plastic)]() [![GitHub forks](https://img.shields.io/github/forks/stfbnc/fathon.svg?style=social&label=Fork&style=plastic)]() [![Python 3.5+](https://img.shields.io/badge/python-3.5+-blue.svg)](https://www.python.org/) [![Generic badge](https://img.shields.io/badge/status-active-pink.svg)](https://shields.io/) [![DOI](https://zenodo.org/badge/214290119.svg)](https://zenodo.org/badge/latestdoi/214290119) [![DOI](https://joss.theoj.org/papers/10.21105/joss.01828/status.svg)](https://doi.org/10.21105/joss.01828)

*fathon* is a python package for DFA (*Detrended Fluctuation Analysis*) and related algorithms.

It is mostly written in Cython and C in order to speed up computations.

## Getting started

### Prerequisites

:bangbang: At the moment, *fathon* is only available for Linux and macOS. :bangbang:

To make *fathon* work, you need:

1. **_Python version:_**
  
   - Python 3.5 or higher
   
2. **_Python packages:_**
  - wget
  
   - numpy (version >= 1.15)
   - cython
  
3. **_OpenMP:_**
   - Linux
     - a C compiler with OpenMP support is required (gcc recommended)
   - macOS
     - the default C compiler (via Xcode) does not support OpenMP and *fathon* will generally run slower than on Linux. If you want to exploit parallelisation, you have to install another C compiler that supports OpenMP (like gcc). You also have to modifiy the file <code>setup.py</code> present in the *fathon* package:
       - see [this](https://stackoverflow.com/questions/54776301/cython-prange-is-repeating-not-parallelizing) and links therein for how to set the C compiler
       - inside <code>setup.py</code>, in function <code>get_extension()</code>, replace the code after the <code>if current_os == "Darwin":</code> with the code after the <code>elif current_os == "Linux":</code>

4. **_External libraries:_**

   - The [GSL (*GNU Scientific Library*)](https://www.gnu.org/software/gsl/) is required for *fathon* to work. It is used to perform polynomial fits during the fluctuation function evaluation in the external C functions. It is not mandatory to install GSL by yourself since it will be downloaded and installed in your <code>/usr/local/</code> directory during *fathon* installation if you do not have it (see **Installing** section for more informations).

### Installing

*fathon* can be installed via <code>pip</code>.

*fathon* needs GSL and the paths to GSL includes and libraries to work.

To install *fathon*, you can choose one of the following (**all paths must be absolute**):

- If you do not have GSL installed, first install it by yourself (from source code or via package manager) in your preferred location (say <code>/my/gsl/location/</code>), then run <code>GSLINC="/my/gsl/location/include/" GSLLIB="/my/gsl/location/lib/" pip install fathon</code> to install *fathon*;
- If you do not have GSL installed, run <code>pip install fathon</code>. It will download and install GSL under the <code>/usr/local/</code> directory and then it will install *fathon*. These procedure could take a few minutes;
- If you already have GSL installed, run <code>GSLINC="/path/to/include/" GSLLIB="/path/to/lib/" pip install fathon</code> to install *fathon*. In case you do not know exactly where GSL is installed, you can run <code>find / -name gsl</code> to find the location of <code>/path/to/include/</code> and <code>/path/to/lib/</code>.

:warning: GSL stores the <code>.h</code> files to be included under the <code>.../include/gsl/</code> directory, but the <code>GSLINC</code> path must not go beyond <code>.../include/</code>!

:warning::warning: If you are using a virtual environment you need to add the GSL libraries to your <code>LD_LIBRARY_PATH</code>:

1. Activate the virtual environment;
2. Run <code>LD_LIBRARY_PATH="/your/GSLLIB/path/:$LD_LIBRARY_PATH"</code>. If you let pip install GSL (the second of the installation options previously listed) <code>/your/GSLLIB/path/</code> is <code>/usr/local/lib/</code>.
3. Run <code>export LD_LIBRARY_PATH</code>;
4. Install *fathon*.

## Examples

A jupyter notebook can be used (<code>fathon_example.ipynb</code>) to play with the four algorithms provided by the *fathon* package, namely

- <code>DFA</code> (*Detrended Fluctuation Analysis*)
- <code>MFDFA</code> (*Multifractal Detrended Fluctuation Analysis*)
- <code>DCCA</code> (*Detrended Cross-Correlation Analysis*)
- <code>HT</code> (*Time-dependent Hurst exponent*)

<code>MFDFA</code> also provides methods for mass exponent τ and multifractal spectrum *f*(α) computation.

<code>DCCA</code> has methods for computing the cross-correlation coefficient ρ_DCCA and the corresponding confidence intervals.

If you cannot use the notebook, four Python scripts are provided, <code>dfa.py</code>, <code>mfdfa.py</code>, <code>dcca.py</code>, and <code>ht.py</code>.

Algorithms are implemented on two time series of gaussian white noise, but you can replace them with any time series you like.

## Documentation

[API documentation](https://fathon.readthedocs.io/)

## Contributing

To report bugs or improvements, or for any other question, please see [CONTRIBUTING.md](https://github.com/stfbnc/fathon/blob/master/CONTRIBUTING.md).

##Credits

If you are using *fathon* in your research, please cite:

Bianchi, S., (2020). fathon: A Python package for a fast computation of  detrendend fluctuation analysis and related algorithms. Journal of Open  Source Software, 5(45), 1828, https://doi.org/10.21105/joss.01828

## Version

fathon v0.1.2

## Author

- Stefano Bianchi
  - github - [stfbnc](https://github.com/stfbnc)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](https://github.com/stfbnc/fathon/blob/master/LICENSE) file for details.

This code uses GSL which is licensed under the GNU General Public License v3.0, and can be obtained at https://www.gnu.org/software/gsl/.