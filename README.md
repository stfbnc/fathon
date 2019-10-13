# fathon

*fathon* is a python package for DFA (*Detrended Fluctuation Analysis*) and related algorithms.

It is mostly written in Cython and C in order to speed up computations.

## Getting started

To download the code you have two options:

1. From Terminal, run <code>git clone https://github.com/stfbnc/fathon.git</code> 
2. [Go to the repository](https://github.com/stfbnc/fathon.git) and download the code

### Prerequisites

To make *fathon* work, you need

- Python 3
- GSL (*GNU Scientific Library*), to perform polynomial fits. You can download GSL [here](https://www.gnu.org/software/gsl/)

Since it is possible to exploit parallel computation, a C compiler with OpenMP support is required under Linux (gcc recommended).

Under macOS, the default C compiler does not support OpenMP. To exploit parallelisation, you have to install (even though not recommended) a C compiler that supports OpenMP, like gcc. See for instance [this](https://stackoverflow.com/questions/54776301/cython-prange-is-repeating-not-parallelizing) and links therein.

### Installing

1. Open Terminal and go to the folder where you have downloaded *fathon*
2. Enter the *fathon* folder
3. Run <code>python  setup.py build_ext -- =/path/to/GSL/includes -- =/path/to/GSL/libs --inplace</code>, where <code>/path/to/GSL/includes</code> and <code>/path/to/GSL/libs</code> are the paths where you have installed the GSL includes and libraries respectively
4. In the *fathon* folder there are now four files corresponding to the compiled Cyhon code
5. Move the folder *fathon* to the standard location of Python packages (e.g. <code>/your/path/pythonversion/site-packages</code>) or to where you want to use the package
6. You can now use *fathon* by typing <code>import fathon</code> from Python

## Tests

A jupyter notebook can be used (<code>fathon_example.ipynb</code>) to test the four algorithms provided by the *fathon* package, namely

- <code>DFA</code> (*Detrended Fluctuation Analysis*)
- <code>MFDFA</code> (*Multifractal Detrended Fluctuation Analysis*)
- <code>DCCA</code> (*Detrended Cross-Correlation Analysis*)
- <code>HT</code> (*Time-dependent Hurst exponent*)

<code>MFDFA</code> also provides methods for mass exponent τ and multifractal spectrum *f*(α) computation.

<code>DCCA</code> has methods for computing the cross-correlation coefficient ρ_DCCA and the corresponding confidence intervals.

If you cannot use the notebook, four Python scripts are provided, <code>dfa.py</code>, <code>mfdfa.py</code>, <code>dcca.py</code>, and <code>ht.py</code>.

Algorithms are tested on two time series of gaussian white noise, but you can replace them with any time series you like.

## Version

fathon v0.1

## Support

For bugs or improvements, write to fathon.package@gmail.com

## Author

- Stefano Bianchi
  - github - [stfbnc](https://github.com/stfbnc)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](https://github.com/stfbnc/fathon/blob/master/LICENSE) file for details.
