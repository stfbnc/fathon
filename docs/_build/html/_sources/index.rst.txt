.. fathon documentation master file, created by
   sphinx-quickstart on Wed Nov 20 11:37:07 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

fathon (v0.1.2)
***************

Python package for detrended fluctuation analysis (DFA) and related algorithms.

Current version is available for Linux and macOS only.

Requirements
============

- Python 3.5+
- wget
- numpy (>=1.15)
- cython
- C compiler:

  - On Linux, a C compiler with OpenMP support is required (gcc recommended);
  - On macOS, the default C compiler does not support OpenMP and *fathon* does not exploit parallelisation. Some instructions on possible ways to exploit parallelisation on macOS are given in the README of the source code repository (https://github.com/stfbnc/fathon.git).

 - GSL (https://www.gnu.org/software/gsl/) is not mandatory but nice to have before *fathon* installation (see **Installation** section for further details).

Installation
============

*fathon* can be installed via :code:`pip`.

*fathon* needs GSL and the paths to GSL includes and libraries to work.

To install *fathon*, you can choose one of the following (**all paths must be absolute**):

- If you do not have GSL installed, first install it by yourself (from source code or via package manager) in your preferred location (say :code:`/my/gsl/location/`), then run :code:`GSLINC="/my/gsl/location/include/" GSLLIB="/my/gsl/location/lib/" pip install fathon` to install *fathon*;
- If you do not have GSL installed, run :code:`pip install fathon`. It will download and install GSL under the :code:`/usr/local/` directory and then it will install *fathon*. These procedure could take a few minutes;
- If you already have GSL installed, run :code:`GSLINC="/path/to/include/" GSLLIB="/path/to/lib/" pip install fathon` to install *fathon*. In case you do not know exactly where GSL is installed, you can run :code:`find / -name gsl` to find the location of :code:`/path/to/include/` and :code:`/path/to/lib/`.

**Warning!!!** GSL stores the :code:`.h` files to be included under the :code:`.../include/gsl/` directory, but the :code:`GSLINC` path must not go beyond :code:`.../include/`

**Warning!!!** If you are using a virtual environment you need to add the GSL libraries to your :code:`LD_LIBRARY_PATH`:

- Activate the virtual environment;
- Run :code:`LD_LIBRARY_PATH="/your/GSLLIB/path/:$LD_LIBRARY_PATH"`. If you let pip install GSL (the second of the installation options previously listed) :code:`/your/GSLLIB/path/` is :code:`/usr/local/lib/`.
- Run :code:`export LD_LIBRARY_PATH`;
- Install *fathon*.

Documentation for the Code
==========================
.. toctree::
   :maxdepth: 1

   fun_class/fathon.subtractMean
   fun_class/fathon.toAggregated
   fun_class/fathon.DFA
   fun_class/fathon.MFDFA
   fun_class/fathon.DCCA
   fun_class/fathon.HT
