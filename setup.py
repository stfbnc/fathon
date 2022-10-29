from setuptools import find_packages
from distutils.core import setup, Extension
import os
from Cython.Build import cythonize
import numpy
import platform
import sys
import re


if platform.system() == "Darwin":
    if platform.processor() == "arm":
        os.environ["CC"] = "/usr/local/opt/llvm/bin/clang"
        os.environ["LDFLAGS"] = "-L/usr/local/opt/llvm/lib"
        os.environ["CPPFLAGS"] = "-I/usr/local/opt/llvm/include"
    else:
        os.environ["CC"] = "gcc-11"
        os.environ["CXX"] = "g++-11"

    gsl_inc = "/usr/local/include"
    gsl_lib = "/usr/local/lib/"
elif platform.system() == "Linux":
    pass
elif platform.system() == "Windows":
    gsl_inc = "C:\\gsl\\include"
    gsl_lib = "C:\\gsl\\lib"
else:
    raise ValueError("Cannot build wheel for OS: {}".format(platform.system()))


def get_extension(module_name, src_name, current_os):
    sources = [src_name, os.path.join("fathon", "cLoops.c")]

    if current_os == "Darwin":
        return Extension(module_name,
                         sources=sources,
                         include_dirs=[numpy.get_include(), gsl_inc],
                         library_dirs=[gsl_lib],
                         libraries=["gsl", "gslcblas", "m"],
                         extra_compile_args=["-O2", "-fopenmp"],
                         extra_link_args=["-fopenmp"])
    elif current_os == "Linux":
        return Extension(module_name,
                         sources=sources,
                         include_dirs=[numpy.get_include()],
                         libraries=["gsl", "gslcblas", "m"],
                         extra_compile_args=["-O2", "-fopenmp"],
                         extra_link_args=["-fopenmp"])
    elif current_os == "Windows":
        return Extension(module_name,
                         sources=sources,
                         include_dirs=[numpy.get_include(), gsl_inc],
                         library_dirs=[gsl_lib],
                         libraries=["gsl", "gslcblas"],
                         extra_compile_args=["/O2", "-openmp"],
                         extra_link_args=[])


if __name__ == "__main__":
    if sys.version_info[0] == 3:
        running_os = platform.system()
        
        extensions = [get_extension("fathon.dfa", os.path.join("fathon", "dfa.pyx"), running_os),
                      get_extension("fathon.dcca", os.path.join("fathon", "dcca.pyx"), running_os),
                      get_extension("fathon.mfdfa", os.path.join("fathon", "mfdfa.pyx"), running_os),
                      get_extension("fathon.mfdcca", os.path.join("fathon", "mfdcca.pyx"), running_os),
                      get_extension("fathon.ht", os.path.join("fathon", "ht.pyx"), running_os)]
                          
        README = ""
        chk = 0
        readme_file = open(os.path.join("docs", "index.rst"), "r", encoding="utf8")
        for line in readme_file:
            if line[:6] == "fathon":
                chk = 1
            if line[:13] == "Documentation":
                chk = 0
            if chk == 1:
                README += re.sub(":code:", "", line)
        readme_file.close()

        setup(name="fathon",
              version="1.3.2",
              author="Stefano Bianchi",
              author_email="fathon.package@gmail.com",
              url="https://github.com/stfbnc/fathon.git",
              license="GPLv3.0",
              description="A pyhton package for detrended fluctuation analysis (DFA) and related algorithms.",
              long_description_content_type="text/markdown",
              long_description=README,
              packages=find_packages(),
              classifiers=["License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
                           "Operating System :: MacOS",
                           "Operating System :: Unix",
                           "Operating System :: Microsoft :: Windows",
                           "Programming Language :: Cython",
                           "Programming Language :: C",
                           "Programming Language :: Python :: 3.7",
                           "Programming Language :: Python :: 3.8",
                           "Programming Language :: Python :: 3.9",
                           "Programming Language :: Python :: 3.10",
                           "Topic :: Scientific/Engineering"],
              python_requires=">=3.7",
              install_requires=["numpy>=1.20"],
              project_urls={"Documentation": "https://fathon.readthedocs.io/",
                            "Bug Reports": "https://github.com/stfbnc/fathon/issues",
                            "Source": "https://github.com/stfbnc/fathon/"},
              ext_modules=cythonize(extensions),
              package_data={"fathon": ["LICENSE"]},
              include_package_data=True)
    else:
        sys.exit("fathon requires python 3.")
