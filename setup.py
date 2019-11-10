from setuptools import find_packages
from distutils.core import setup, Extension
import subprocess
import os
from Cython.Build import cythonize
import numpy
import platform
import sys
from pathlib import Path

fathon_path = sys.path[-1]+"/fathon"
gsl_path = fathon_path+"/fathonGSL"

def gsl_install():
    for p in sys.path:
        process = subprocess.Popen("if [ -d "+p+"/fathon ]; then rm -rf "+p+"/fathon; fi", shell=True, cwd=Path(__file__).parent.absolute())
        process.wait()
    command = "mkdir "+fathon_path+" && mkdir "+gsl_path+" && cd src/gsl_code/ && ./configure --prefix="+gsl_path+" && make && make install && cd .. && rm -rf gsl_code"
    process = subprocess.Popen(command, shell=True, cwd=Path(__file__).parent.absolute())
    process.wait()

def get_extension(module_name, src_name, current_os):
    sources = [src_name, "src/cLoops.c"]
    include_dirs = [numpy.get_include(), gsl_path+"/include"]
    library_dirs = [gsl_path+"/lib"]
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

def move_fathon():
    mv_fathon = "cp __init__.py "+fathon_path+" && cp tsHelper.py "+fathon_path+" && cp README.md "+fathon_path+" && cp LICENSE "+fathon_path+" && mv dcca* "+fathon_path+" && mv dfa* "+fathon_path+" && mv ht* "+fathon_path+" && mv mfdfa* "+fathon_path
    process = subprocess.Popen(mv_fathon, shell=True, cwd=Path(__file__).parent.absolute())
    process.wait()

if __name__ == '__main__':
    if sys.version_info[0] == 3:
        running_os = platform.system()
        if running_os != "Windows":
            gsl_install()
        
            extensions = [get_extension("dfa", "src/dfa.pyx", running_os),
                          get_extension("dcca", "src/dcca.pyx", running_os),
                          get_extension("mfdfa", "src/mfdfa.pyx", running_os),
                          get_extension("ht", "src/ht.pyx", running_os)]

            setup(name="fathon",
                  version="0.1",
                  author="Stefano Bianchi",
                  url="https://github.com/stfbnc/fathon.git",
                  license="GPLv3.0",
                  description="pyhton package for detrended fluctuation analysis (DFA) and related algorithms.",
                  packages=find_packages(),
                  install_requires=["numpy", "cython"],
                  ext_modules=cythonize(extensions, build_dir="build")
                  )
        
            move_fathon()
        else:
            print("fathon is not available on Windows yet.")
    else:
        sys.exit("fathon requires pyhton 3.")
