"""
CSV-driven bin classes for HISP reactor modeling.
Each bin represents one row from the CSV configuration table.

Naming conventions
------------------
- **flux_id** (formerly *Bin number*): identifies the physical location / bin
  in the binned‑flux‑data file (``Bin_Index`` column).  Many CSV rows can
  share the same ``flux_id`` (e.g. sub‑bins of the same parent geometry bin).
- **sim_id** (formerly *bin_id*): a **unique** identifier for each simulation
  row.  If the CSV contains a *Sim. ID* column the value is taken from there;
  otherwise the 1‑based row number is used automatically.
"""

from typing import Optional, Any
from dataclasses import dataclass
import pandas as pd
from materials.materials import Material


@dataclass
class BinConfiguration:
    """Configuration parameters for HISP simulation."""
    rtol: float
    atol: float
    fp_max_stepsize: float  # FP max. stepsize (s)
    max_stepsize_no_fp: float  # Max. stepsize no FP (s)
    bc_plasma_facing_surface: str  # BC Plasma Facing Surface
    bc_rear_surface: str  # BC rear surface
    
    def __post_init__(self):
        """Validate configuration parameters."""
        if self.rtol <= 0:
            raise ValueError(f"rtol must be positive, got {self.rtol}")
        if self.atol <= 0:
            raise ValueError(f"atol must be positive, got {self.atol}")
        if self.fp_max_stepsize <= 0:
            raise ValueError(f"fp_max_stepsize must be positive, got {self.fp_max_stepsize}")
        if self.max_stepsize_no_fp <= 0:
            raise ValueError(f"max_stepsize_no_fp must be positive, got {self.max_stepsize_no_fp}")


class Bin:
    """
    A bin class representing one row from the CSV configuration table.
    Each bin contains all geometric, material, and simulation properties.
    """
    
    def __init__(
        self,
        flux_id: int,
        material: Material,
        thickness: float,
        cu_thickness: float,
        mode: str,
        parent_bin_surf_area: float,
        surface_area: float,
        f_ion_flux_fraction: float,
        location: str,
        z_start: float = 0.0,
        r_start: float = 0.0,
        z_end: float = 0.0,
        r_end: float = 0.0,
        coolant_temp: float = 343.0,
        bin_configuration: Optional[BinConfiguration] = None,
        sim_id: Optional[int] = None,
        calculate_implantation_params: bool = True,
        atom_view_factor: float = 1.0,
        # ── legacy aliases (ignored if the new names are provided) ──
        bin_number: Optional[int] = None,
        bin_id: Optional[int] = None,
    ):
        """
        Initialize a CSV-based bin.
        
        Args:
            flux_id: Flux ID – links this row to the ``Bin_Index`` column in
                the binned-flux-data files.  (Formerly called *Bin number*.)
            z_start: Z coordinate start position (m)  [optional]
            r_start: R coordinate start position (m)  [optional]
            z_end: Z coordinate end position (m)  [optional]
            r_end: R coordinate end position (m)  [optional]
            material: Material object (W, B, SS, etc.)
            thickness: Bin thickness (m)
            cu_thickness: Copper thickness (m)
            mode: Operating mode (hw, lw, shadowed, wetted, etc.)
            parent_bin_surf_area: Surface area of parent bin (m^2)
            surface_area: Surface area of this specific bin/mode (m^2)
            f_ion_flux_fraction: Ion flux fraction
            location: Location identifier (FW, DIV, etc.)
            coolant_temp: Coolant temperature (K)
            bin_configuration: BinConfiguration object with simulation parameters
            sim_id: Unique simulation identifier.  If *None*, defaults to
                ``flux_id``.
            atom_view_factor: Scalar multiplier applied to the **atomic** flux
                for this bin (default 1.0 = no change).
            calculate_implantation_params: whether to compute implantation
                parameters from flux data at runtime.
            
        Calculated Properties:
            ion_scaling_factor: f_ion_flux_fraction * parent_bin_surf_area / surface_area
        """
        # ── handle legacy keyword aliases ──
        if bin_number is not None and flux_id is None:
            flux_id = bin_number
        if bin_id is not None and sim_id is None:
            sim_id = bin_id

        # Geometric properties (optional – may be 0 when not provided)
        self.z_start = z_start
        self.r_start = r_start
        self.z_end = z_end
        self.r_end = r_end
        
        # Material properties
        if not isinstance(material, Material):
            raise TypeError(
                "Bin.material must be a Material instance. "
                "Ensure CSVBinLoader matches material names to materials.csv"
            )
        self.material = material
        self.thickness = thickness
        self.cu_thickness = cu_thickness
        
        # Operating properties
        self.mode = mode
        self.parent_bin_surf_area = parent_bin_surf_area
        self.surface_area = surface_area
        self.f_ion_flux_fraction = f_ion_flux_fraction
        self.location = location
        self.coolant_temp = coolant_temp
        
        # ── identifiers ──
        self.flux_id = flux_id
        self.sim_id = sim_id if sim_id is not None else flux_id

        # ── backward-compatible aliases ──
        # Kept so that any old code doing ``bin.bin_number`` or ``bin.bin_id``
        # still works without changes.
        self.bin_number = self.flux_id
        self.bin_id = self.sim_id

        # Atom view factor
        self.atom_view_factor = atom_view_factor

        # Calculated properties
        self.ion_scaling_factor = self.f_ion_flux_fraction * self.parent_bin_surf_area / self.surface_area

        # Simulation configuration (use provided or create default)
        self.bin_configuration = bin_configuration if bin_configuration is not None else BinConfiguration(
            rtol=1e-10,
            atol=1e10,
            fp_max_stepsize=5.0,
            max_stepsize_no_fp=100.0,
            bc_plasma_facing_surface="Robin - Surf. Rec. + Implantation",
            bc_rear_surface="Neumann - no flux"
        )
        
        # Implantation parameters (computed at runtime from plasma data)
        self.implantation_params = None
        
        # Control flag for calculating implantation parameters from flux data
        self.calculate_implantation_params = calculate_implantation_params

    @property
    def material_name(self) -> str:
        """Return the material name string (for legacy code expecting a string)."""
        return getattr(self.material, 'name', str(self.material))
    
    @property
    def copper_thickness(self) -> float:
        """Compatibility property for HISP temperature functions."""
        return self.cu_thickness
    
    @property
    def start_point(self) -> tuple[float, float]:
        """Get start point as (Z, R) tuple."""
        return (self.z_start, self.r_start)
    
    @property
    def end_point(self) -> tuple[float, float]:
        """Get end point as (Z, R) tuple."""
        return (self.z_end, self.r_end)
    
    @property
    def length(self) -> float:
        """Calculate the poloidal length of the bin (m)."""
        return (
            (self.z_end - self.z_start) ** 2 + 
            (self.r_end - self.r_start) ** 2
        ) ** 0.5
    
    @property
    def is_first_wall(self) -> bool:
        """Check if this is a first wall bin."""
        return self.location.upper() == "FW"
    
    @property
    def is_divertor(self) -> bool:
        """Check if this is a divertor bin."""
        return self.location.upper() in ["DIV", "DIVERTOR"]
    
    @property
    def is_shadowed(self) -> bool:
        """Check if this bin is in shadowed mode."""
        return self.mode.lower() in ["shadowed", "shadow"]
    
    @property
    def is_wetted(self) -> bool:
        """Check if this bin is in any wetted mode."""
        return self.mode.lower() in ["wetted", "wet", "hw", "lw", "high_wetted", "low_wetted"]
    
    def get_implantation_data(self, pulse, plasma_data_handling, ion: bool = True):
        """
        Extract implantation data (energy, angle) for a specific pulse and particle type.
        
        Args:
            pulse: Pulse object from scenario
            plasma_data_handling: PlasmaDataHandling object with flux data
            ion: True for ions, False for atoms
            
        Returns:
            Dictionary with 'energy' and 'angle' keys (values may be None if not available)
        """
        energy = None
        angle = None
        
        try:
            pulse_type = getattr(pulse, 'pulse_type', None)
            if pulse_type is None:
                return {'energy': energy, 'angle': angle}
            
            pulse_data = plasma_data_handling.pulse_type_to_data.get(pulse_type)
            if pulse_data is None:
                return {'energy': energy, 'angle': angle}
            
            # flux_id matches Bin_Index in the flux data files
            bin_row = pulse_data[pulse_data['Bin_Index'] == self.flux_id]
            
            if bin_row.empty:
                return {'energy': energy, 'angle': angle}
            
            if ion:
                energy = bin_row['E_ion'].values[0]
                angle = bin_row['alpha_ion'].values[0]
            else:
                energy = bin_row['E_atom'].values[0]
                angle = bin_row['alpha_atom'].values[0]
            
            if pd.isna(energy):
                energy = None
            if pd.isna(angle):
                angle = None
                
        except (KeyError, IndexError, AttributeError):
            pass
        
        return {
            'energy': energy,
            'angle': angle
        }
    
    def __str__(self) -> str:
        """String representation of the bin."""
        return (
            f"Bin(sim_id={self.sim_id}, flux_id={self.flux_id}, "
            f"material={self.material_name}, mode={self.mode}, "
            f"location={self.location}, thickness={self.thickness*1000:.1f}mm)"
        )
    
    def __repr__(self) -> str:
        """Detailed representation of the bin."""
        return self.__str__()


class BinCollection:
    """Collection of CSV-based bins."""
    
    def __init__(self, bins: list[Bin] = None):
        """Initialize collection with list of Bin objects."""
        self.bins = bins if bins is not None else []
    
    def add_bin(self, bin: Bin):
        """Add a bin to the collection."""
        self.bins.append(bin)
    
    def get_bin_by_sim_id(self, sim_id: int) -> Bin:
        """Get bin by its simulation ID."""
        for bin in self.bins:
            if bin.sim_id == sim_id:
                return bin
        raise ValueError(f"No bin found with sim_id {sim_id}")

    # Backward-compatible alias
    def get_bin_by_id(self, sim_id: int) -> Bin:
        """Alias for get_bin_by_sim_id (backward compatibility)."""
        return self.get_bin_by_sim_id(sim_id)
    
    def get_bin_by_flux_id(self, flux_id: int) -> Bin:
        """Get first bin matching the given flux_id."""
        for bin in self.bins:
            if bin.flux_id == flux_id:
                return bin
        raise ValueError(f"No bin found with flux_id {flux_id}")

    # Backward-compatible alias
    def get_bin_by_number(self, flux_id: int) -> Bin:
        """Alias for get_bin_by_flux_id (backward compatibility)."""
        return self.get_bin_by_flux_id(flux_id)
    
    def get_bins_by_material(self, material: str) -> list[Bin]:
        """Get all bins with specified material."""
        return [bin for bin in self.bins if bin.material_name.upper() == material.upper()]
    
    def get_bins_by_location(self, location: str) -> list[Bin]:
        """Get all bins at specified location (FW, DIV, etc.)."""
        return [bin for bin in self.bins if bin.location.upper() == location.upper()]
    
    def get_bins_by_mode(self, mode: str) -> list[Bin]:
        """Get all bins with specified mode."""
        return [bin for bin in self.bins if bin.mode.lower() == mode.lower()]
    
    @property
    def first_wall_bins(self) -> list[Bin]:
        """Get all first wall bins."""
        return [bin for bin in self.bins if bin.is_first_wall]
    
    @property
    def divertor_bins(self) -> list[Bin]:
        """Get all divertor bins."""
        return [bin for bin in self.bins if bin.is_divertor]
    
    def __len__(self) -> int:
        """Return number of bins in collection."""
        return len(self.bins)
    
    def __iter__(self):
        """Make collection iterable."""
        return iter(self.bins)
    
    def __str__(self) -> str:
        """String representation of the collection."""
        fw_count = len(self.first_wall_bins)
        div_count = len(self.divertor_bins)
        return f"BinCollection({len(self.bins)} bins: {fw_count} FW, {div_count} DIV)"


class Reactor(BinCollection):
    """
    A reactor representing the complete collection of all bins from a CSV table.
    This is the main class for representing the entire ITER reactor configuration.
    """
    
    def __init__(self, bins: list[Bin] = None, csv_path: str = None):
        """
        Initialize reactor with list of CSVBin objects.
        
        Args:
            bins: List of Bin objects representing all reactor bins
            csv_path: Optional path to the source CSV file for reference
        """
        super().__init__(bins)
        self.csv_path = csv_path
    
    @classmethod
    def from_csv(cls, csv_path: str) -> 'Reactor':
        """
        Create a complete reactor by loading all bins from a CSV table.
        
        Args:
            csv_path: Path to the CSV configuration file
            
        Returns:
            Reactor: Complete reactor with all bins from the CSV table
        """
        from bins_from_csv.csv_bin_loader import CSVBinLoader
        
        loader = CSVBinLoader(csv_path)
        bin_collection = loader.load_all_bins()
        return cls(bins=bin_collection.bins, csv_path=csv_path)
    
    @property
    def total_bins(self) -> int:
        """Get total number of bins in the reactor."""
        return len(self.bins)
    
    @property
    def materials_summary(self) -> dict[str, int]:
        """Get summary of materials used in the reactor."""
        materials = {}
        for bin in self.bins:
            material = bin.material_name.upper()
            materials[material] = materials.get(material, 0) + 1
        return materials
    
    @property
    def locations_summary(self) -> dict[str, int]:
        """Get summary of bin locations in the reactor."""
        locations = {}
        for bin in self.bins:
            location = bin.location.upper()
            locations[location] = locations.get(location, 0) + 1
        return locations
    
    def get_reactor_summary(self) -> str:
        """Get comprehensive summary of the reactor configuration."""
        summary = [
            f"Reactor Summary:",
            f"  Total bins: {self.total_bins}",
            f"  CSV source: {self.csv_path or 'Not specified'}",
            f"  First Wall bins: {len(self.first_wall_bins)}",
            f"  Divertor bins: {len(self.divertor_bins)}",
            f"  Materials: {self.materials_summary}",
            f"  Locations: {self.locations_summary}"
        ]
        return "\n".join(summary)
    
    def __str__(self) -> str:
        """String representation of the reactor."""
        return f"Reactor({self.total_bins} total bins: {len(self.first_wall_bins)} FW, {len(self.divertor_bins)} DIV)"
