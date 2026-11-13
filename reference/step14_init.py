# REFERENCE for migration Step 14 (not part of REF; untracked). Goes to medicc/fstengine/__init__.py.
# REF's medicc/native/__init__.py with native -> fstengine. Same as REF otherwise.

"""fstengine: the OpenFST copy-number engine (C++ in src/, ctypes wrapper in model.py).

P = A + ADA + ADADA with A = L_LOH o G (losses, including LOH, followed by gains)
and the capped whole-genome doubling D(c) = min(2c, 8). Every event costs 1, copy
numbers are capped at 8 and ancestors are reconstructed without pruning. Allele-specific
and total copy-number data use the same model.
"""
from medicc.fstengine.model import Model, ModelConfig
from medicc.fstengine.scheduler import pairwise_distances, reconstruct
