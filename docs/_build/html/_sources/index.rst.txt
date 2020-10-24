.. fathon documentation master file, created by
   sphinx-quickstart on Wed Nov 20 11:37:07 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

fathon (v1.1)
***************

Current version is available for Linux (64 bit) and macOS only.

Python package for detrended fluctuation analysis (DFA) and related algorithms.
`fathon` provides four main algorithms, namely

- `DFA` (Detrended Fluctuation Analysis)
- `MFDFA` (Multifractal Detrended Fluctuation Analysis)
- `DCCA` (Detrended Cross-Correlation Analysis)
- `HT` (Time-dependent Hurst exponent)

`MFDFA` also provides methods for the mass exponent τ and the multifractal spectrum f(α).

`DCCA` has methods for computing the cross-correlation coefficient ρ_DCCA and the corresponding confidence intervals.

Requirements
============

- Python 3.5+
- numpy (>=1.15)
- Cython

Installation
============

`pip install fathon`

Documentation for the Code
==========================
.. toctree::
   :maxdepth: 1

   fun_class/fathon.fathonUtils
   fun_class/fathon.DFA
   fun_class/fathon.MFDFA
   fun_class/fathon.DCCA
   fun_class/fathon.HT
