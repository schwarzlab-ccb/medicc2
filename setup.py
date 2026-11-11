# Package metadata lives in pyproject.toml. This file only defines the Cython/C++ extensions.
# medicc.native._engine is a plain C++ shared library (loaded with ctypes, not imported).
import glob
import os
import sys

from Cython.Build import cythonize
from setuptools import Extension, setup

include_dirs = ['fstlib/cext']
library_dirs = []
extra_link_args = []
conda_path = os.environ.get("CONDA_PREFIX", None)
if conda_path:
  include_dirs.append(os.path.join(conda_path, 'include'))
  library_dirs.append(os.path.join(conda_path, 'lib'))
  extra_link_args.append('-Wl,-rpath,{}'.format(os.path.join(conda_path, 'lib')))

extra_compile_args = ['-std=c++17']
if sys.platform.startswith("darwin"):
  extra_compile_args.append("-stdlib=libc++")
  extra_compile_args.append("-mmacosx-version-min=10.14")

fst_libraries = ["fst", "fstfar", "fstscript", "fstfarscript"]

setup(
    ext_modules=cythonize([
        Extension("fstlib.cext.pywrapfst",
                  ["fstlib/cext/pywrapfst.pyx"],
                  include_dirs=include_dirs,
                  library_dirs=library_dirs,
                  libraries=fst_libraries,
                  extra_compile_args=extra_compile_args,
                  extra_link_args=extra_link_args,
                  language="c++"),

        Extension("fstlib.cext.ops",
                  ["fstlib/cext/ops.pyx"],
                  include_dirs=include_dirs,
                  library_dirs=library_dirs,
                  libraries=fst_libraries,
                  extra_compile_args=extra_compile_args,
                  extra_link_args=extra_link_args,
                  language="c++")
    ], include_path=['fstlib/cext']) + [
        Extension("medicc.native._engine",
                  ["medicc/native/src/native.cpp"],
                  depends=glob.glob("medicc/native/src/*.hpp"),
                  include_dirs=include_dirs,
                  library_dirs=library_dirs,
                  libraries=["fst"],
                  extra_compile_args=extra_compile_args + ['-O3'],
                  extra_link_args=extra_link_args,
                  language="c++")
    ]
)
