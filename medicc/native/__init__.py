"""Native OpenFST copy-number model.

P = A + ADA + ADADA with A = L_LOH o G (losses, including LOH, followed by gains)
and the capped whole-genome doubling D(c) = min(2c, 8). Every event costs 1, copy
numbers are capped at 8 and ancestors are reconstructed without pruning. Allele-specific
and total copy-number data use the same model (see README).
"""
from medicc.native.model import Model, ModelConfig
from medicc.native.scheduler import pairwise_distances, reconstruct
