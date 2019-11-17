from setuptools import find_packages
from distutils.core import setup, Extension
import subprocess
import os
from Cython.Build import cythonize
import numpy
import platform
import sys
import urllib.request

main_path = os.path.dirname(os.path.abspath(__file__))
home_path = os.path.expanduser("~")

def gsl_install():
	process = subprocess.Popen(os.path.join(".", "fathon", "fathon_gsl_install"), cwd=main_path)
	process.wait()
	if process.returncode != 0:
		sys.exit("Failed to install GSL.")

gsl_inc = os.environ.get("GSLINC", None)
gsl_lib = os.environ.get("GSLLIB", None)
if gsl_inc is None and gsl_lib is None:
    urllib.request.urlretrieve("ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz", os.path.join(home_path, "gsl-latest.tar.gz"))
    gsl_install()
    gsl_inc = "/usr/local/include"
    gsl_lib = "/usr/local/lib"
elif gsl_inc is not None and gsl_lib is not None:
	pass
else:
	sys.exit("Both GSLINC and GSLLIB must or must not be given.")

def get_extension(module_name, src_name, current_os):
    sources = [src_name, os.path.join("fathon", "cLoops.c")]
    include_dirs = [numpy.get_include(), gsl_inc]
    library_dirs = [gsl_lib]
    libraries = ["gsl", "gslcblas", "m"]
    extra_compile_args_macos = ["-O2"]
    extra_compile_args_linux = ["-O2", "-fopenmp"]
    extra_link_args = ["-fopenmp"]
    if current_os == "Darwin":
        return Extension(module_name,
                         sources=sources,
                         include_dirs=include_dirs,
                         library_dirs=library_dirs,
                         libraries=libraries,
                         extra_compile_args=extra_compile_args_macos)
    elif current_os == "Linux":
        return Extension(module_name,
                         sources=sources,
                         include_dirs=include_dirs,
                         library_dirs=library_dirs,
                         libraries=libraries,
                         extra_compile_args=extra_compile_args_linux,
                         extra_link_args=extra_link_args)

if __name__ == "__main__":
    if sys.version_info[0] == 3:
        running_os = platform.system()
        if running_os != "Windows":
            extensions = [get_extension("fathon.dfa", os.path.join("fathon", "dfa.pyx"), running_os),
                          get_extension("fathon.dcca", os.path.join("fathon", "dcca.pyx"), running_os),
                          get_extension("fathon.mfdfa", os.path.join("fathon", "mfdfa.pyx"), running_os),
                          get_extension("fathon.ht", os.path.join("fathon", "ht.pyx"), running_os)]

            setup(name="fathon",
                  version="0.1.2",
                  author="Stefano Bianchi",
                  author_email="fathon.package@gmail.com",
                  url="https://github.com/stfbnc/fathon.git",
                  license="GPLv3.0",
                  description="pyhton package for detrended fluctuation analysis (DFA) and related algorithms.",
                  packages=find_packages(),
                  install_requires=["numpy", "cython"],
                  ext_modules=cythonize(extensions),
                  package_data={"fathon": ["LICENSE"]},
                  include_package_data=True)
        else:
            sys.exit("fathon is not available on Windows yet.")
    else:
        sys.exit("fathon requires python 3.")
