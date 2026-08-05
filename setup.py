import sys
import os
from pathlib import Path

from Cython.Build import cythonize
from setuptools import Extension, setup

sys.path.append('fstlib/cext')
conda_path = os.environ.get("CONDA_PREFIX", None)
extra_link_args = ['-L{}/lib'.format(conda_path)]

extra_compile_args = ['-std=c++17']
if sys.platform.startswith("darwin"):
  extra_compile_args.append("-stdlib=libc++")
  extra_compile_args.append("-mmacosx-version-min=10.14")

# Define __version__
__version__ = None
exec(open("medicc/_version.py").read())

setup(
    name='medicc2',
    # Change this in medicc/_version.py
    version=__version__,
    author='Tom L Kaufmann, Marina Petkovic, Chenxi Nie, Roland F Schwarz',
    author_email='tkau93@gmail.com, marina.55kovic@gmail.com, roland.f.schwarz@gmail.com',
    description='Whole-genome doubling-aware copy number phylogenies for cancer evolution',
    long_description=(Path(__file__).parent / "README.md").read_text(),
    long_description_content_type='text/markdown',
    url='https://bitbucket.org/schwarzlab/medicc2',
    classifiers=[
                "Programming Language :: Python :: 3",
                "License :: OSI Approved :: MIT License",
                "Operating System :: OS Independent",
    ],
    packages=['medicc', 'fstlib', 'fstlib.cext'],
    scripts=['medicc2'],
    license='GPL-3',
    install_requires=[
        'numpy>=1.20.1',
        'pyyaml>=5.4.1',
        'pandas>=1.2',
        'biopython>=1.78',
        'scipy>=1.7',
        'matplotlib>=3.3.4'
    ],
    extras_require={
        'Parallelization': ['joblib>=1.0.1'],
    },
    package_data={
        "medicc": ["objects/*.fst", "objects/*.bed", "logging_conf.yaml"],
        "fstlib": ["logging_conf.yaml", "cext/*.pxd", "cext/*.pyx", "cext/*.h"],
    },
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
                  language = "c++"),
        Extension(
            "medicc.acyclic_prune",
            ["medicc/cext/acyclic_prune.pyx"],
            include_dirs=["medicc/cext", "fstlib/cext"],
            depends=["medicc/cext/acyclic_prune.h"],
            libraries=["fst", "fstfar", "fstscript", "fstfarscript"],
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args,
            language="c++",
        )
    ])
)
