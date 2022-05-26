# fathon ![Linux](https://github.com/stfbnc/fathon/actions/workflows/linux.yml/badge.svg) ![macOS](https://github.com/stfbnc/fathon/actions/workflows/macos.yml/badge.svg) ![Windows](https://ci.appveyor.com/api/projects/status/tl2a8c84bbvxu37p/branch/master?svg=true&passingText=Windows&pendingText=Windows&failingText=Windows)

[![Issues](https://img.shields.io/github/issues-raw/stfbnc/fathon.svg?maxAge=25000)](https://github.com/stfbnc/fathon/issues) [![GitHub stars](https://img.shields.io/github/stars/stfbnc/fathon.svg?style=social&label=Stars&style=plastic)]() [![GitHub forks](https://img.shields.io/github/forks/stfbnc/fathon.svg?style=social&label=Fork&style=plastic)]() [![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/)

[![DOI](https://zenodo.org/badge/214290119.svg)](https://zenodo.org/badge/latestdoi/214290119) [![DOI](https://joss.theoj.org/papers/10.21105/joss.01828/status.svg)](https://doi.org/10.21105/joss.01828)



`fathon` is a python package for DFA (*Detrended Fluctuation Analysis*) and related algorithms.

It is mostly written in Cython and C in order to speed up computations.

`fathon` provides five main algorithms, namely

- <code>DFA</code> (*Detrended Fluctuation Analysis*)
- <code>MFDFA</code> (*Multifractal Detrended Fluctuation Analysis*)
- <code>DCCA</code> (*Detrended Cross-Correlation Analysis*)
- <code>MFDCCA</code> (*Multifractal Detrended Cross-Correlation Analysis*)
- <code>HT</code> (*Time-dependent Hurst exponent*)

<code>MFDFA</code> and <code>MFDCCA</code> also provides methods for the mass exponent τ and the multifractal spectrum *f*(α).

<code>DCCA</code> has methods to compute the cross-correlation coefficient ρ_DCCA and the corresponding confidence intervals.

### Supported platforms

|    Linux x86_64    |    Linux ARM64     |    macOS x86_64    | macOS ARM64 |   Windows 64bit    |
|:------------------:| :----------------: | :----------------: |:-----------:| :----------------: |
| :white_check_mark: | :white_check_mark: | :white_check_mark: |     :x:     | :white_check_mark: |

### Prerequisites

 - Python 3.7 or higher
 - numpy (version >= 1.20)

### Installing

As easy as `pip install fathon`

### Examples

A jupyter notebook can be used (<code>fathon_example.ipynb</code>) to play with the five algorithms provided by the `fathon` package.

If you cannot use the notebook, five Python scripts are provided, <code>dfa.py</code>, <code>mfdfa.py</code>, <code>dcca.py</code>, <code>mfdcca.py</code>, and <code>ht.py</code>.

Algorithms are implemented on two time series of gaussian white noise, but you can replace them with any time series you like.

## Documentation [![Documentation Status](https://readthedocs.org/projects/fathon/badge/?version=latest)](https://fathon.readthedocs.io/en/latest/?badge=latest)

[API documentation](https://fathon.readthedocs.io/)

## Contributing

To report bugs or improvements, or for any other question, please see [CONTRIBUTING.md](https://github.com/stfbnc/fathon/blob/master/CONTRIBUTING.md).

## Credits

If you are using `fathon` in your research, please cite:

Bianchi, S., (2020). fathon: A Python package for a fast computation of  detrendend fluctuation analysis and related algorithms. Journal of Open  Source Software, 5(45), 1828, https://doi.org/10.21105/joss.01828

## Version  [![PyPI version](https://badge.fury.io/py/fathon.svg)](https://badge.fury.io/py/fathon)

fathon v1.3.1

## Changelog

#### v1.3.1

- faster algorithms

#### v1.3

- <code>MFDCCA</code> algorithm
- <code>overlap</code> option for <code>DCCA</code>, to allow using both overlapping and non-overlapping windows
- OpenMP also for Windows

#### v1.2

- few adjustments to C extensions for Windows' C compiler compatibility

#### v1.1

- save object state to binary file and reload it later

#### v1.0

- wheels! :ferris_wheel::ferris_wheel:
- no more pre-installing step of the GSL library :tada::tada:
- window's sizes array must be now passed to all the methods
- `logBase` option for methods that perform fits
- `verbose` option
- pre-computed `hq0` can be now passed to the `computeHT` method

#### v0.1.2

- first release

## Author

- Stefano Bianchi
  - github - [stfbnc](https://github.com/stfbnc)

## License  [![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-orange.svg)](https://opensource.org/licenses/)

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](https://github.com/stfbnc/fathon/blob/master/LICENSE) file for details.

This code uses GSL which is licensed under the GNU General Public License v3.0, and can be obtained at https://www.gnu.org/software/gsl/.
