from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy
import platform

def get_extension(module_name, src_name, current_os):
    sources = [src_name, "src/cLoops.c"]
    include_dirs = [numpy.get_include()]
    libraries = ["gsl", "gslcblas", "m"]
    extra_compile_args_macos = ["-O2"]
    extra_compile_args_linux = ["-O2", "-fopenmp"]
    extra_link_args = ["-fopenmp"]
    if current_os == "Darwin":
        return Extension(module_name,
                         sources=sources,
                         include_dirs=include_dirs,
                         libraries=libraries,
                         extra_compile_args=extra_compile_args_macos)
    elif current_os == "Linux":
        return Extension(module_name,
                         sources=sources,
                         include_dirs=include_dirs,
                         libraries=libraries,
                         extra_compile_args=extra_compile_args_linux,
                         extra_link_args=extra_link_args)

running_os = platform.system()
extensions = [get_extension("dfa", "src/dfa.pyx", running_os),
              get_extension("dcca", "src/dcca.pyx", running_os),
              get_extension("mfdfa", "src/mfdfa.pyx", running_os),
              get_extension("ht", "src/ht.pyx", running_os)]

setup(
      name="fathon",
      ext_modules=cythonize(extensions, build_dir="build")
)
