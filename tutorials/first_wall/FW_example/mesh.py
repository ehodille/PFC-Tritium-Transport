"""
mesh.py — Tutorial: First Wall
==============================
Defines BINS_MESHES: a graded mesh for each bin, concentrated at the
plasma-facing surface where implantation occurs.

See wiki/Input-Files-Reference.md for mesh configuration details.
"""

import os
import sys
import numpy as np
from pathlib import Path
from typing import Dict

# ── Ensure repo root is on sys.path ──────────────────────────────────────────
# When run by the cluster runner, PFC_TT_PATH is already set.
# When developing locally, add the repo root manually if needed.
_repo_root = Path(__file__).resolve().parents[3]  # tutorials/first_wall/FW_example/ → repo root
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

from meshing import MeshBin
from bins_from_csv.csv_bin_loader import CSVBinLoader


def graded_vertices(L: float, h0: float, r: float) -> np.ndarray:
    """
    Generate graded mesh vertices from 0 to L.

    Parameters
    ----------
    L  : domain length (thickness) in metres
    h0 : first element size at x=0 (plasma-facing surface)
    r  : growth ratio > 1  (elements grow toward x=L)

    Returns
    -------
    1-D numpy array of vertex positions
    """
    xs = [0.0]
    h = h0
    while xs[-1] + h < L:
        xs.append(xs[-1] + h)
        h *= r
    if xs[-1] < L:
        xs.append(float(L))
    return np.array(xs, dtype=float)


def symmetric_graded_vertices(L: float, h0: float, r: float) -> np.ndarray:
    """
    Generate a symmetric graded mesh refined at both x=0 and x=L.

    Grading grows from each surface toward the centre; the two halves
    are merged and deduplicated.  Use this when both surfaces carry an
    active boundary condition (e.g. Dirichlet at plasma face + Dirichlet
    at rear for permeation studies).

    Parameters
    ----------
    L  : domain length (thickness) in metres
    h0 : first element size at each surface
    r  : growth ratio > 1

    Returns
    -------
    1-D numpy array of vertex positions
    """
    half = L / 2.0
    # left half: grow from x=0 toward centre
    left = [0.0]
    h = h0
    while left[-1] + h < half:
        left.append(left[-1] + h)
        h *= r
    # right half: mirror of left
    right = [L - x for x in left]
    right.reverse()
    pts = np.unique(np.array(left + right + [L], dtype=float))
    return pts


# ── Load reactor to access bin thicknesses ────────────────────────────────────
# INPUT_DIR_CONTEXT is set by the SLURM runner to the absolute input folder.
# Falls back to the directory containing this file when run locally.
_input_dir = os.environ.get("INPUT_DIR_CONTEXT", str(Path(__file__).parent))
_csv_path = os.path.join(_input_dir, "input_table.csv")
_mat_path = os.path.join(_input_dir, "materials.csv")

_loader = CSVBinLoader(_csv_path, materials_csv_path=_mat_path)
_reactor = _loader.load_reactor()

# ── Build BINS_MESHES ─────────────────────────────────────────────────────────
BINS_MESHES: Dict[int, MeshBin] = {}

for _bin in _reactor.bins:
    _bc_plasma = _bin.bin_configuration.bc_plasma_facing_surface

    if "Dirichlet" in _bc_plasma:
        # Dirichlet analytical implantation approximation at the plasma-facing
        # surface: concentration prescribed at x=0.  Use the standard one-sided
        # graded mesh refined at the surface.
        _mesh = graded_vertices(L=_bin.thickness, h0=1e-8, r=1.02)
    else:
        # Robin (surface recombination + implantation) or similar: flux-based BC
        # couples both surfaces to the interior.  Use a symmetric mesh refined
        # at both x=0 and x=L so neither end is under-resolved.
        _mesh = symmetric_graded_vertices(L=_bin.thickness, h0=1e-10, r=1.01)

    BINS_MESHES[_bin.sim_id] = MeshBin(sim_id=_bin.sim_id, mesh=_mesh)

print(f"✓ Tutorial mesh: generated meshes for {len(BINS_MESHES)} bin(s)")
