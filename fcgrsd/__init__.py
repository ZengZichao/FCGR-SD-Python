"""FCGR-SD: alignment-free genome comparison on Grassmann manifolds.

Python implementation of

    Li X., Zhou T., Feng X., Yau S.-T., Yau S. S.-T. (2024).
    Exploring geometry of genome space via Grassmann manifolds.
    The Innovation 5(5): 100677.  https://doi.org/10.1016/j.xinn.2024.100677

Modules
-------
``cgr``        Chaos Game Representation coordinates and grid indices
``fcgr``       Frequency matrix of Chaos Game Representation (FCGR)
``grassmann``  Truncated SVD, principal angles, generalized Grassmann distance
``distance``   Pairwise distance matrices
``tree``       Neighbor-joining and Newick I/O
``metrics``    Robinson-Foulds distance and classical MDS
``io``         FASTA and distance-matrix I/O
``pipeline``   End-to-end workflow
``cli``        Command-line entry point
"""

from __future__ import annotations

from .cgr import ALPHABETS, cgr_coordinates, cgr_grid_indices
from .distance import distance_matrix
from .fcgr import fcgr_matrix
from .grassmann import (
    DEFAULT_EPS,
    grassmann_distance,
    principal_angles,
    subspace_basis,
    truncated_basis,
)
from .metrics import classical_mds, robinson_foulds
from .tree import Tree, neighbor_joining

__version__ = "1.0.0"

__all__ = [
    "ALPHABETS",
    "DEFAULT_EPS",
    "Tree",
    "__version__",
    "cgr_coordinates",
    "cgr_grid_indices",
    "classical_mds",
    "distance_matrix",
    "fcgr_matrix",
    "grassmann_distance",
    "neighbor_joining",
    "principal_angles",
    "robinson_foulds",
    "subspace_basis",
    "truncated_basis",
]
