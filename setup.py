"""Build a standalone MEDICC branch, including local OpenFST and full-FST engines."""
from pathlib import Path
import os
import shutil
import subprocess
import sys
from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext
from Cython.Build import cythonize
ROOT = Path(__file__).resolve().parent
PREFIX = Path(os.getenv('CONDA_PREFIX', sys.prefix))
include = [str(ROOT), str(ROOT/'fstlib/cext'), str(ROOT/'medicc/cext'), str(PREFIX/'include')]
extensions = [Extension(name, [source], include_dirs=include,
    library_dirs=[str(PREFIX/'lib')], runtime_library_dirs=[str(PREFIX/'lib')],
    libraries=['fst', 'fstfar', 'fstscript', 'fstfarscript'], extra_compile_args=['-std=c++17'], language='c++')
    for name, source in [('fstlib.cext.pywrapfst', 'fstlib/cext/pywrapfst.pyx'),
                         ('fstlib.cext.ops', 'fstlib/cext/ops.pyx'),
                         ('medicc.acyclic_prune', 'medicc/cext/acyclic_prune.pyx')]]

class BuildEngines(build_ext):
    def run(self):
        super().run()
        destination = ROOT if self.inplace else Path(self.build_lib)
        for name, source, library in [('ancestors','ancestors.cpp','libdirect_wgd2.so'),
                                     ('direct_pair','pairwise.cpp','libpair_direct.so'),
                                     ('reference_pair','pairwise.cpp','libpair_direct.so')]:
            folder = ROOT/'medicc/cext/full_fst'/name
            target = destination/'medicc/cext/full_fst'/name/'build'
            target.mkdir(parents=True, exist_ok=True)
            shutil.copy2(folder/'assets/ancestral_parts.txt', target/'ancestral_parts.txt')
            command = [os.getenv('CXX','g++'), '-O3','-std=c++17','-pthread','-shared','-fPIC',
                '-DANCESTRAL_MAX_WGDS=2',str(folder/'src'/source)]
            if name == 'ancestors':
                command += [str(ROOT/'medicc/cext/full_fst/common/src/openfst_language.cpp'),
                            f'-I{PREFIX}/include', f'-L{PREFIX}/lib', f'-Wl,-rpath,{PREFIX}/lib', '-lfst']
            subprocess.run(command + ['-o',str(target/library)], check=True)

setup(name='medicc2-full-fst', version='1.1.2.dev2',
      description='MEDICC reference implementation of full L_LOH with up to two WGDs',
      license='GPL-3.0', packages=find_packages(), python_requires='>=3.10',
      install_requires=['numpy','pandas','scipy','biopython','matplotlib','pyyaml','joblib'],
      extras_require={'events': ['pyranges>=0.1,<0.2']},
      entry_points={'console_scripts':['medicc2=medicc.full_fst.cli:main']},
      package_data={'medicc':[
                        'logging_conf.yaml', 'objects/*.bed', 'cext/*.h', 'cext/*.pyx',
                        'cext/full_fst/*/assets/*.txt', 'cext/full_fst/*/src/*.cpp',
                        'cext/full_fst/*/src/*.hpp', 'cext/full_fst/*/src/generated/*.hpp',
                        'cext/full_fst/*/build/*.so', 'cext/full_fst/*/build/*.txt'],
                    'fstlib':['logging_conf.yaml','cext/*.h','cext/*.pxd','cext/*.pyx']},
      ext_modules=cythonize(extensions, include_path=include, build_dir='build/cython', language_level=3),
      cmdclass={'build_ext':BuildEngines})
