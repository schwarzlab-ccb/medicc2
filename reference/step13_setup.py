# REFERENCE for migration Step 13 (not part of REF; untracked). Goes to setup.py in NEW.
# Step 12's setup.py with the metadata removed: it now lives in pyproject.toml ([project]).
# setup.py only defines the compiled extensions (fstlib's Cython modules and the fstengine).
import glob  # Step 12: the header list for `depends`
import sys
import os

from Cython.Build import cythonize
from setuptools import Extension, setup

sys.path.append('fstlib/cext')
conda_path = os.environ.get("CONDA_PREFIX", None)
extra_link_args = ['-L{}/lib'.format(conda_path)]

extra_compile_args = ['-std=c++17']
if sys.platform.startswith("darwin"):
  extra_compile_args.append("-stdlib=libc++")
  extra_compile_args.append("-mmacosx-version-min=10.14")

setup(
    ext_modules = cythonize([
        Extension("fstlib.cext.pywrapfst", 
                  ["fstlib/cext/pywrapfst.pyx"],
                  include_dirs=['fstlib/cext'],
                  libraries=["fst", "fstfar", "fstscript", "fstfarscript"],
                  extra_compile_args=extra_compile_args,
                  extra_link_args=extra_link_args,
                  language = "c++"),

        Extension("fstlib.cext.ops", 
                  ["fstlib/cext/ops.pyx"],
                  include_dirs=['fstlib/cext'],
                  libraries=["fst", "fstfar", "fstscript", "fstfarscript"],
                  extra_compile_args=extra_compile_args,
                  extra_link_args=extra_link_args,
                  language = "c++")
    ]) + [
        # Step 12: the fstengine. Plain C++, so it stays outside cythonize([...]).
        # One translation unit: fstengine.cpp #includes every header in src/.
        # `depends`: rebuild when any header changes, not only fstengine.cpp.
        # Builds to medicc/fstengine/_fstengine.cpython-312-x86_64-linux-gnu.so (model.py finds it).
        Extension("medicc.fstengine._fstengine",
                  ["medicc/fstengine/src/fstengine.cpp"],
                  depends=glob.glob("medicc/fstengine/src/*.hpp"),
                  libraries=["fst"],
                  extra_compile_args=extra_compile_args + ['-O3'],
                  extra_link_args=extra_link_args,
                  language="c++"),
    ]
)
