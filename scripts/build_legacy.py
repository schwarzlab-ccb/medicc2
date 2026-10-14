"""Local build of the copied wrappers; does not install or edit any sibling package."""
from pathlib import Path
import os
import sys
from setuptools import Extension, setup
from Cython.Build import cythonize
ROOT = Path(__file__).resolve().parents[1]
PREFIX = Path(os.getenv('CONDA_PREFIX', sys.prefix))
FST = ROOT/'fstlib/cext'
MED = ROOT/'medicc/cext'
include = [str(ROOT), str(FST), str(MED), str(PREFIX/'include')]
extensions = [Extension(name, [str(source)], include_dirs=include,
    library_dirs=[str(PREFIX/'lib')], runtime_library_dirs=[str(PREFIX/'lib')],
    libraries=['fst', 'fstfar', 'fstscript', 'fstfarscript'],
    extra_compile_args=['-std=c++17'], language='c++')
    for name, source in [('fstlib.cext.pywrapfst', FST/'pywrapfst.pyx'),
                         ('fstlib.cext.ops', FST/'ops.pyx'),
                         ('medicc.acyclic_prune', MED/'acyclic_prune.pyx')]]
setup(name='medicc-full-fst-local-extensions', version='0.1',
      ext_modules=cythonize(extensions, include_path=include, build_dir=str(ROOT/'build/cython'), language_level=3))
