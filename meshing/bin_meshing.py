import numpy as np
from typing import Dict


class MeshBin:
    """Container linking a sim_id to its mesh array."""
    def __init__(self, sim_id: int = None, mesh: np.ndarray = None, *, bin_id: int = None):
        """
        Initialize with sim_id and mesh array.
        
        Args:
            sim_id: Simulation ID of the bin
            mesh: Array of mesh vertices
            bin_id: **Deprecated** alias for *sim_id* (backward compat)
        """
        if sim_id is None and bin_id is not None:
            sim_id = bin_id
        if sim_id is None:
            raise ValueError("MeshBin requires a sim_id (or legacy bin_id)")
        self.sim_id = sim_id
        # backward-compatible alias
        self.bin_id = sim_id
        self.mesh = mesh
    
    def __repr__(self):
        return f"MeshBin(sim_id={self.sim_id}, nodes={len(self.mesh)})"

