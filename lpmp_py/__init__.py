from .raw_solvers import gm_solver, mgm_solver

try:
    import torch
except ImportError as e:
    import warnings
    warnings.warn("Missing Torch, skipping import of torch modules")
else:
    from .torch_wrappers import *
