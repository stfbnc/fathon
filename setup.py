from setuptools import find_packages
from distutils.core import setup, Extension
#import subprocess
import os
from Cython.Build import cythonize
import numpy
import platform
import sys
import re

main_path = os.path.dirname(os.path.abspath(__file__))
home_path = os.path.expanduser("~")

#def gsl_install():
#	process = subprocess.Popen(os.path.join(".", "fathon", "fathon_gsl_install"), cwd=main_path)
#	process.wait()
#	if process.returncode != 0:
#		sys.exit("Failed to install GSL.")

#gsl_inc = os.environ.get("GSLINC", None)
#gsl_lib = os.environ.get("GSLLIB", None)
#if gsl_inc is None and gsl_lib is None:
#    wget.download("ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz", os.path.join(home_path, "gsl-latest.tar.gz"))
#gsl_install()


iGSL_1 = "chmod 777 fathon_gsl_install"
os.system(iGSL_1)
iGSL_2 = "./fathon_gsl_install"
os.system(iGSL_2)


gsl_inc = "./fathon/3rd_party/gsl/include"
gsl_lib = "./fathon/3rd_party/gsl/lib/"
omp_inc = "./fathon/3rd_party/omp/include"
omp_lib = "./fathon/3rd_party/omp/lib/"
#    gsl_inc = "/usr/local/include"
#    gsl_lib = "/usr/local/lib"

#if running_os == "Darwin":
#    cmd1 = "export CFLAGS=\"-Xpreprocessor -fopenmp $CFLAGS\""
#    cmd2 = "export CXXFLAGS=\"-Xpreprocessor -fopenmp $CXXFLAGS\""
#    os.system(cmd1)
#    os.system(cmd2)

def get_extension(module_name, src_name, current_os):
    sources = [src_name, os.path.join("fathon", "cLoops.c")]

    if current_os == "Darwin":
        cmd1 = "install_name_tool -id \"@loader_path/3rd_party/gsl/lib/libgslcblas.dylib\" " + gsl_lib + "libgslcblas.dylib"
        cmd2 = "install_name_tool -id \"@loader_path/3rd_party/gsl/lib/libgsl.dylib\" " + gsl_lib + "libgsl.dylib"
        os.system(cmd1)
        os.system(cmd2)
        return Extension(module_name,
                         sources=sources,
                         include_dirs=[numpy.get_include(), gsl_inc],#, omp_inc],
                         library_dirs=[gsl_lib],#, omp_lib],
                         libraries=["gsl", "gslcblas", "m"],
                         extra_compile_args=["-O2", "-fopenmp"],
                         extra_link_args=["-fopenmp"])
                         #extra_link_args=["-lomp"])
                         #runtime_library_dirs=["@rpath/3rd_party/gsl/lib/"],
                         #extra_objects=[gsl_lib+"libgsl.a", gsl_lib+"libgslcblas.a"])
    elif current_os == "Linux":
        return Extension(module_name,
                         sources=sources,
                         include_dirs=[numpy.get_include(), gsl_inc],
                         library_dirs=[gsl_lib],
                         runtime_library_dirs=["$ORIGIN/3rd_party/gsl/lib/"],
                         libraries=["gsl", "gslcblas", "m"],
                         extra_compile_args=["-O2", "-fopenmp"],
                         extra_link_args=["-fopenmp"])

if __name__ == "__main__":
    if sys.version_info[0] == 3:
        running_os = platform.system()
        if running_os != "Windows":
            extensions = [get_extension("fathon.dfa", os.path.join("fathon", "dfa.pyx"), running_os),
                          get_extension("fathon.dcca", os.path.join("fathon", "dcca.pyx"), running_os),
                          get_extension("fathon.mfdfa", os.path.join("fathon", "mfdfa.pyx"), running_os),
                          get_extension("fathon.ht", os.path.join("fathon", "ht.pyx"), running_os)]
                          
            README = ""
            chk = 0
            readme_file = open("docs/index.rst", "r")
            for line in readme_file:
                if line[:6] == "fathon":
                    chk = 1
                if line[:13] == "Documentation":
                    chk = 0
                if chk == 1:
                    README += re.sub(":code:", "", line)
            readme_file.close()

            setup(name="fathon",
                  version="0.1.2.post4",
                  author="Stefano Bianchi",
                  author_email="fathon.package@gmail.com",
                  url="https://github.com/stfbnc/fathon.git",
                  license="GPLv3.0",
                  description="pyhton package for detrended fluctuation analysis (DFA) and related algorithms.",
                  long_description_content_type="text/markdown",
                  long_description=README,
                  packages=find_packages(),
                  classifiers=["License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
                               "Operating System :: MacOS",
                               "Operating System :: Unix",
                               "Programming Language :: Cython",
                               "Programming Language :: C",
                               "Programming Language :: Python :: 3.5",
                               "Programming Language :: Python :: 3.6",
                               "Programming Language :: Python :: 3.7",
                               "Programming Language :: Python :: 3.8",
                               "Topic :: Scientific/Engineering"],
                  python_requires=">=3.5",
                  install_requires=["numpy>=1.15", "Cython"],
                  project_urls={"Documentation": "https://fathon.readthedocs.io/",
                                "Bug Reports": "https://github.com/stfbnc/fathon/issues",
                                "Source": "https://github.com/stfbnc/fathon/"},
                  ext_modules=cythonize(extensions),
                  package_data={"fathon": ["LICENSE"]},
                  include_package_data=True)
        else:
            sys.exit("fathon is not available on Windows yet.")
    else:
        sys.exit("fathon requires python 3.")
