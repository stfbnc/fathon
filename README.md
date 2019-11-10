# fathon

*fathon* is a python package for DFA (*Detrended Fluctuation Analysis*) and related algorithms.

It is mostly written in Cython and C in order to speed up computations.

## Getting started

### Prerequisites

At the moment, *fathon* is only available for Linux and macOS.

To make *fathon* work, you need

- Python 3
- Cython

on both Linux and macOs.

Under Linux only, a C compiler with OpenMP support is required (gcc recommended).

Under macOS, the default C compiler (via Xcode) does not support OpenMP and *fathon* will generally run slower than on Linux. If you want to exploit parallelisation, you have to install another C compiler that supports OpenMP (like gcc). You also have to modifiy the file <code>setup.py</code> present in the *fathon* package:

- see [this](https://stackoverflow.com/questions/54776301/cython-prange-is-repeating-not-parallelizing) and links therein for how to set the C compiler
- inside <code>setup.py</code>, in function <code>get_extension()</code>, replace the code after the <code>if current_os == "Darwin":</code> with the code after the <code>elif current_os == "Linux":</code>

### Installing

1. Open the Terminal, go to the directory where you want to download the *fathon* package, and run <code>git clone https://github.com/stfbnc/fathon.git</code>
2. <code>cd fathon</code>
3. Run <code>python  setup.py build_ext --inplace</code>. This will take a few minutes. It will
   - create a package named <code>fathon</code> in the last path of your current Python's <code>sys.path</code>
   - install the [GSL (*GNU Scientific Library*)](https://www.gnu.org/software/gsl/) in a folder named <code>fathonGSL</code> inside the package <code>fathon</code>
   - compile the Cython scripts
   - move all the necessary files to the package <code>fathon</code>
4. You can now use *fathon* by typing <code>import fathon</code> from Python

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

## Version

fathon v0.1

## Support

For bugs or improvements, open an issue at https://github.com/stfbnc/fathon.git

## Author

- Stefano Bianchi
  - github - [stfbnc](https://github.com/stfbnc)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](https://github.com/stfbnc/fathon/blob/master/LICENSE) file for details.

This code uses GSL which is licensed under the GNU General Public License v3.0, and can be obtained at https://www.gnu.org/software/gsl/.