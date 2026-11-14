# REFERENCE for migration Step 13 (not part of REF's code). Goes to setup.py in NEW.
# Step 12's setup.py plus two changes, marked "Step 13": medicc.fstengine in `packages`, and its
# files in `package_data`. pyproject.toml stays as it is: the metadata already lives here.
# (REF moved all metadata into pyproject.toml; that is a clean-up, not needed to build or ship.)
import glob  # Step 12: the header list for `depends`
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
    packages=['medicc', 'medicc.fstengine', 'fstlib', 'fstlib.cext'],  # Step 13: + medicc.fstengine
    scripts=['medicc2'],
    license='GPL-3',
    python_requires='>=3.12',
    install_requires=[
        'numpy>=2.0,<3',
        'pyyaml>=5.4.1',
        'pandas>=3.0,<4',
        'biopython>=1.78',
        'scipy>=1.7',
        'matplotlib>=3.3.4'
    ],
    extras_require={
        'Parallelization': ['joblib>=1.0.1'],
    },
    package_data={
        # Step 13: + objects/*.txt (ancestral_parts.txt, read by the engine at run time)
        "medicc": ["objects/*.fst", "objects/*.bed", "objects/*.txt", "logging_conf.yaml"],
        # Step 13: the engine sources, so a source install can build the extension
        "medicc.fstengine": ["src/*.hpp", "src/*.cpp"],
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
