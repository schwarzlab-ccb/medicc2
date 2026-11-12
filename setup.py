# Package metadata lives in pyproject.toml. This file only defines the Cython/C++ extensions.
# medicc.native._engine is a plain C++ shared library (loaded with ctypes, not imported).
import glob
import os
import subprocess
import sys
from pathlib import Path

from Cython.Build import cythonize
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

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

# GCC profile-guided optimization applies only to the native engine.
# Both passes must use the same source, compiler and build directory.
native_compile_args = extra_compile_args + ['-O3']
native_link_args = list(extra_link_args)
pgo_mode = os.environ.get('MEDICC_PGO', '')
if pgo_mode not in ('', 'generate', 'use'):
  raise RuntimeError("MEDICC_PGO must be 'generate' or 'use', or unset")
if pgo_mode:
  if sys.platform == 'darwin':
    # Reserved for the macOS/Clang profile generation and merge workflow.
    raise RuntimeError("PGO for macOS is not implemented yet; unset MEDICC_PGO "
                       "to use the standard build")
  elif not sys.platform.startswith('linux'):
    raise RuntimeError(f"PGO is not implemented for {sys.platform}; unset "
                       "MEDICC_PGO to use the standard build")

  # Linux/GCC implementation.
  pgo_dir = Path(os.environ.get('MEDICC_PGO_DIR', 'build/pgo')).resolve()
  if pgo_mode == 'generate':
    pgo_dir.mkdir(parents=True, exist_ok=True)
    pgo_args = [f'-fprofile-generate={pgo_dir}', '-fprofile-update=atomic']
  else:
    if not any(pgo_dir.rglob('*.gcda')):
      raise RuntimeError(f"No GCC profiles found in {pgo_dir}; build with "
                         "MEDICC_PGO=generate and run a training workload first")
    pgo_args = [f'-fprofile-use={pgo_dir}', '-Werror=missing-profile',
                '-Werror=coverage-mismatch']
  native_compile_args += pgo_args
  native_link_args += pgo_args


class NativeBuildExt(build_ext):
  """Check the configured C++ compiler before a PGO build."""

  def build_extensions(self):
    if pgo_mode:
      compiler = getattr(self.compiler, 'compiler_cxx', None)
      if not compiler:
        raise RuntimeError("PGO requires a GCC C++ compiler; set CC=gcc and CXX=g++")
      try:
        result = subprocess.run(compiler + ['-dM', '-E', '-x', 'c++', '-'],
                                input='', text=True, capture_output=True,
                                check=True, timeout=30)
      except (OSError, subprocess.SubprocessError) as exc:
        raise RuntimeError("Could not identify the C++ compiler for PGO; "
                           "set CC=gcc and CXX=g++") from exc
      macros = {line.split()[1] for line in result.stdout.splitlines()
                if line.startswith('#define ') and len(line.split()) >= 2}
      other_compilers = {'__clang__', '__INTEL_COMPILER',
                         '__INTEL_LLVM_COMPILER', '__NVCOMPILER'}
      if '__GNUC__' not in macros or macros & other_compilers:
        raise RuntimeError("PGO currently supports GCC on Linux only; "
                           "set CC=gcc and CXX=g++ or unset MEDICC_PGO")
    super().build_extensions()


setup(
    cmdclass={'build_ext': NativeBuildExt},
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
                  extra_compile_args=native_compile_args,
                  extra_link_args=native_link_args,
                  language="c++")
    ]
)
