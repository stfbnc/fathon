from setuptools.command.install import install
from distutils.core import setup, Extension
from setuptools import find_packages
import subprocess
import os
from Cython.Build import cythonize
import numpy
import platform
from pathlib import Path

home_path = str(Path.home()) ##3.5+ altrimenti from os.path import expanduser - home = expanduser("~")

def gsl_install():
    command = "mkdir -p "+home_path+"/fathonGSL && cd src/gsl_code/ && ./configure --prefix="+home_path+"/fathonGSL && make && make install && cd .. && rm -rf gsl_code"
    process = subprocess.Popen(command, shell=True, cwd=Path(__file__).parent.absolute())
    process.wait()

#class CustomInstall(install):
#    def run(self):
#        command = "mkdir -p "+home_path+"/fathonGSL && cd src/gsl_code/ && ./configure --prefix="+home_path+"/fathonGSL && make && make install && cd .. && rm -rf gsl_code"
#        process = subprocess.Popen(command, shell=True, cwd=Path(__file__).parent.absolute())
#        process.wait()
#        install.run(self)

def get_extension(module_name, src_name, current_os):
    sources = [src_name, "src/cLoops.c"]
    include_dirs = [numpy.get_include(), home_path+"/fathonGSL/include"]
    library_dirs = [home_path+"/fathonGSL/lib"]
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

if __name__ == '__main__':
    gsl_install()

    running_os = platform.system()
    if running_os != "Windows":
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
    else:
        print("fathon is not available on Windows yet.")

#cmdclass={'install': CustomInstall},
#include_package_data=True,
