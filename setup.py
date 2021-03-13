from setuptools import find_packages
from distutils.core import setup, Extension
import os
from Cython.Build import cythonize
import numpy
import platform
import sys
import re

if platform.system() == "Darwin":
    os.environ["CC"] = "gcc-9"
    os.environ["CXX"] = "g++-9"

if platform.system() != "Windows":
	iGSL_1 = "chmod 777 fathon_gsl_install"
	os.system(iGSL_1)
	iGSL_2 = "./fathon_gsl_install"
	os.system(iGSL_2)

	gsl_inc = "./fathon/3rd_party/gsl/include"
	gsl_lib = "./fathon/3rd_party/gsl/lib/"
else:
	gsl_inc = "fathon\\3rd_party\\gsl\\include"
	gsl_lib = "fathon\\3rd_party\\gsl\\lib"

def get_extension(module_name, src_name, current_os):
    sources = [src_name, os.path.join("fathon", "cLoops.c")]

    if current_os == "Darwin":
        cmd1 = "install_name_tool -id \"@loader_path/3rd_party/gsl/lib/libgslcblas.dylib\" " + gsl_lib + "libgslcblas.dylib"
        cmd2 = "install_name_tool -id \"@loader_path/3rd_party/gsl/lib/libgsl.dylib\" " + gsl_lib + "libgsl.dylib"
        os.system(cmd1)
        os.system(cmd2)

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
                         include_dirs=[numpy.get_include(), gsl_inc],
                         runtime_library_dirs=["$ORIGIN/3rd_party/gsl/lib/"],
                         libraries=["gsl", "gslcblas", "m"],
                         library_dirs=[gsl_lib],
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
              version="1.3",
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
                           "Topic :: Scientific/Engineering"],
              python_requires=">=3.7",
              install_requires=["numpy>=1.20", "Cython"],
              project_urls={"Documentation": "https://fathon.readthedocs.io/",
                            "Bug Reports": "https://github.com/stfbnc/fathon/issues",
                            "Source": "https://github.com/stfbnc/fathon/"},
              ext_modules=cythonize(extensions),
              package_data={"fathon": ["LICENSE"]},
              include_package_data=True)
    else:
        sys.exit("fathon requires python 3.")
