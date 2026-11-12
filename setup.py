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

# Profile-guided optimization applies only to the native engine.
# Both passes must use the same source, compiler and build directory.
native_compile_args = extra_compile_args + ['-O3']
native_link_args = list(extra_link_args)
pgo_mode = os.environ.get('MEDICC_PGO', '')
if pgo_mode not in ('', 'generate', 'use'):
  raise RuntimeError("MEDICC_PGO must be 'generate' or 'use', or unset")
if pgo_mode:
  if sys.platform != 'darwin' and not sys.platform.startswith('linux'):
    raise RuntimeError(f"PGO is not implemented for {sys.platform}; unset "
                       "MEDICC_PGO to use the standard build")

  pgo_dir = Path(os.environ.get('MEDICC_PGO_DIR', 'build/pgo')).resolve()
  if pgo_mode == 'generate':
    pgo_dir.mkdir(parents=True, exist_ok=True)
    # Clang writes default_%m.profraw (binary-specific, with online merging);
    # GCC writes .gcda files. Atomic counters also support threaded training.
    pgo_args = [f'-fprofile-generate={pgo_dir}', '-fprofile-update=atomic']
  elif sys.platform == 'darwin':
    profile = pgo_dir / 'default.profdata'
    if not profile.is_file():
      raise RuntimeError(f"No merged Clang profile found at {profile}; build with "
                         "MEDICC_PGO=generate, run training workloads, then merge "
                         "the .profraw files with matching llvm-profdata (see README)")
    pgo_args = [f'-fprofile-use={profile}', '-Werror=profile-instr-out-of-date',
                '-Werror=profile-instr-unprofiled']
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
      required = ('Clang C++ compiler; set CC=clang and CXX=clang++'
                  if sys.platform == 'darwin' else
                  'GCC C++ compiler; set CC=gcc and CXX=g++')
      compiler = getattr(self.compiler, 'compiler_cxx', None)
      if not compiler:
        raise RuntimeError(f"PGO requires a {required}")
      try:
        result = subprocess.run(compiler + ['-dM', '-E', '-x', 'c++', '-'],
                                input='', text=True, capture_output=True,
                                check=True, timeout=30)
      except (OSError, subprocess.SubprocessError) as exc:
        raise RuntimeError(f"Could not identify the C++ compiler for PGO; use a {required}") from exc
      macros = {line.split()[1] for line in result.stdout.splitlines()
                if line.startswith('#define ') and len(line.split()) >= 2}
      other_compilers = {'__INTEL_COMPILER', '__INTEL_LLVM_COMPILER', '__NVCOMPILER'}
      if sys.platform == 'darwin':
        if '__clang__' not in macros or macros & other_compilers:
          raise RuntimeError(f"PGO on macOS requires a {required}, or unset MEDICC_PGO")
      elif '__GNUC__' not in macros or macros & (other_compilers | {'__clang__'}):
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
