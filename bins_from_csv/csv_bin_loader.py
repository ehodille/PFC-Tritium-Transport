"""
CSV loader for creating Bin objects from CSV configuration files.

Column matching is **flexible**: spaces, underscores, dots, parentheses,
and letter case are all ignored when looking up column names.  For example,
the following all match the same logical column:

    Flux ID | flux_id | FLUX_ID | flux id | Flux_ID

Each logical column also accepts a set of *aliases* (e.g. the old name
``Bin number`` is still accepted as an alias for ``Flux ID``).
"""

import re
import pandas as pd
from typing import List, Dict, Any, Optional
from pathlib import Path
from bins_from_csv.csv_bin import Bin, BinCollection, Reactor, BinConfiguration
from materials.materials_loader import load_materials


# ---------------------------------------------------------------------------
# Flexible column-name resolution
# ---------------------------------------------------------------------------

def _normalise(name: str) -> str:
    """Collapse a column name to a canonical lowercase form with no
    spaces, underscores, dots, or parentheses."""
    return re.sub(r"[\s_.()\-]+", "", str(name).strip().lower())


# Mapping from *canonical alias* → *internal key* used by the loader.
# Several human-readable forms are listed for every logical column.
_COLUMN_ALIASES: Dict[str, str] = {}

def _register(internal_key: str, *aliases: str):
    """Register one or more human-friendly aliases for an internal key."""
    for alias in aliases:
        _COLUMN_ALIASES[_normalise(alias)] = internal_key

# ── required columns ──
_register("flux_id",
          "Flux ID", "flux_id", "flux id", "FluxID", "Flux_ID",
          "Bin number", "bin_number", "binnumber", "Bin_number")
_register("material",        "Material", "material", "mat")
_register("thickness",       "Thickness (m)", "thickness", "thickness_m", "Thickness")
_register("cu_thickness",    "Cu thickness (m)", "cu_thickness", "cuthickness", "Cu thickness", "Cu_thickness")
_register("mode",            "mode", "Mode")
_register("parent_area",     "S. Area parent bin (m^2)", "s area parent bin", "parent_area",
          "S Area parent bin", "parent area", "sareaparentbin")
_register("surface_area",    "Surface area (m^2)", "surface_area", "surfacearea", "Surface area")
_register("f_ion",           "Ion flux wetted fraction", "ion_flux_wetted_fraction",
          "ionfluxwettedfraction",
          # legacy aliases
          "f (ion flux scaling factor)", "f", "f_ion", "fionfluxscalingfactor",
          "ion flux scaling factor")
_register("location",        "location", "Location", "loc")

# ── optional columns ──
_register("z_start",         "Z_start (m)", "z_start", "zstart", "Z start", "Z_start")
_register("r_start",         "R_start (m)", "r_start", "rstart", "R start", "R_start")
_register("z_end",           "Z_end (m)", "z_end", "zend", "Z end", "Z_end")
_register("r_end",           "R_end (m)", "r_end", "rend", "R end", "R_end")
_register("coolant_temp",    "Coolant Temp. (K)", "coolant_temp", "coolanttemp", "Coolant Temp")
_register("rtol",            "rtol", "Rtol")
_register("atol",            "atol", "Atol")
_register("fp_max_stepsize", "FP max. stepsize (s)", "fp_max_stepsize", "fpmaxstepsize",
          "FP max stepsize", "FP_max_stepsize")
_register("max_stepsize_no_fp", "Max. stepsize no FP (s)", "max_stepsize_no_fp",
          "maxstepsizenofp", "Max stepsize no FP")
_register("bc_plasma",       "BC Plasma Facing Surface", "bc_plasma", "bcplasmafacingsurface",
          "BC Plasma Facing", "bc_plasma_facing_surface")
_register("bc_rear",         "BC rear surface", "bc_rear", "bcrearsurface", "BC_rear_surface")
_register("calc_implant",    "Calculate Implantation Parameters", "calculateimplantationparameters",
          "calc_implant", "Calculate Implantation")
_register("sim_id",          "Sim. ID", "sim_id", "simid", "Sim ID", "SimID",
          "Simulation ID", "simulation_id", "simulationid")
_register("atom_view_factor","Atom view factor", "atom_view_factor", "atomviewfactor",
          "Atom_view_factor", "atom view factor")


def _resolve_column(col_name: str) -> str:
    """Return the internal key for a CSV column, or the original name if unknown."""
    return _COLUMN_ALIASES.get(_normalise(col_name), col_name)


def _build_column_map(df_columns) -> Dict[str, str]:
    """Build a mapping  internal_key → actual_csv_column_name  for every
    column present in the DataFrame."""
    col_map: Dict[str, str] = {}
    for col in df_columns:
        key = _resolve_column(col)
        if key not in col_map:          # first match wins
            col_map[key] = col
    return col_map


class CSVBinLoader:
    """Loads Bin objects from CSV configuration files."""
    
    def __init__(self, csv_path: str, materials_csv_path: Optional[str] = None):
        """
        Initialize loader with CSV file path.
        
        Args:
            csv_path: Path to the CSV configuration file
            materials_csv_path: Optional explicit path to materials CSV
        """
        self.csv_path = csv_path
        # Allow caller to pass either the direct path or the filename located in input_files/
        try:
            self.df = pd.read_csv(csv_path)
        except FileNotFoundError:
            # try under input_files/ for convenience
            alt = Path("input_files") / csv_path
            try:
                self.df = pd.read_csv(alt)
                self.csv_path = str(alt)
            except FileNotFoundError:
                raise FileNotFoundError(f"CSV file not found at '{csv_path}' or '{alt}'")

        # Build flexible column map once
        self._col_map = _build_column_map(self.df.columns)

        self._validate_csv()

        # Load materials CSV if available (default: input_files/materials.csv)
        if materials_csv_path is None:
            materials_csv_path = Path("input_files/materials.csv")
        try:
            self.materials = load_materials(materials_csv_path)
            print(f"✓ Loaded {len(self.materials)} materials from {materials_csv_path}")
        except Exception as e:
            raise RuntimeError(f"Failed to load materials from {materials_csv_path}: {e}")

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------

    def _csv_col(self, internal_key: str) -> Optional[str]:
        """Return the *actual* CSV column name for an internal key, or None."""
        return self._col_map.get(internal_key)

    def _has_col(self, internal_key: str) -> bool:
        return internal_key in self._col_map

    def _validate_csv(self):
        """Validate that CSV has required columns (using flexible matching)."""
        required_keys = [
            "flux_id", "material", "thickness", "cu_thickness", "mode",
            "parent_area", "surface_area", "f_ion", "location",
        ]
        
        missing = [k for k in required_keys if not self._has_col(k)]
        if missing:
            raise ValueError(
                f"Missing required columns in CSV (internal keys): {missing}.\n"
                f"  Recognised columns: {list(self._col_map.keys())}\n"
                f"  CSV header: {list(self.df.columns)}"
            )
        
        print(f"✓ CSV validation passed. Found {len(self.df)} rows with {len(self.df.columns)} columns.")

    def _get(self, row: pd.Series, internal_key: str, default: Any = None) -> Any:
        """Get a value from *row* using flexible column resolution."""
        csv_col = self._csv_col(internal_key)
        if csv_col is None:
            return default
        value = row[csv_col]
        if pd.isna(value):
            return default
        return value

    # kept for any callers that still use the old name
    def _get_column_value(self, row: pd.Series, column: str, default: Any = None) -> Any:
        """Safely get column value – tries flexible resolution first, then literal."""
        key = _resolve_column(column)
        return self._get(row, key, default)

    # ------------------------------------------------------------------
    # bin construction
    # ------------------------------------------------------------------
    
    def load_bin_from_row(self, row: pd.Series, row_index: int) -> Bin:
        """
        Create a Bin from a pandas DataFrame row.
        
        Args:
            row: Pandas Series representing one row from the CSV
            row_index: Index of the row (0-based)
            
        Returns:
            Bin object
        """
        # ── Required columns ──
        flux_id = int(self._get(row, "flux_id"))

        # ── Optional geometry columns (default 0.0) ──
        z_start = float(self._get(row, "z_start", 0.0))
        r_start = float(self._get(row, "r_start", 0.0))
        z_end   = float(self._get(row, "z_end",   0.0))
        r_end   = float(self._get(row, "r_end",   0.0))

        # ── Material ──
        material = str(self._get(row, "material"))
        thickness = float(self._get(row, "thickness"))
        cu_thickness = float(self._get(row, "cu_thickness"))

        mat_obj = None
        mat_name = material.strip()
        if hasattr(self, 'materials') and self.materials:
            mat_obj = self.materials.get(mat_name)
            if mat_obj is None:
                for k, v in self.materials.items():
                    if k.lower() == mat_name.lower():
                        mat_obj = v
                        break
        if mat_obj is None:
            available = ', '.join(sorted(self.materials.keys()))
            raise ValueError(
                f"Unknown material '{material}' in CSV row {row_index + 1}. "
                f"Available materials: {available}"
            )

        # ── Operating properties ──
        mode = str(self._get(row, "mode"))
        parent_bin_surf_area = float(self._get(row, "parent_area"))
        surface_area = float(self._get(row, "surface_area"))
        f_ion_flux_fraction = float(self._get(row, "f_ion"))
        location = str(self._get(row, "location"))

        # ── Optional scalar columns ──
        coolant_temp = float(self._get(row, "coolant_temp", 343.0))

        rtol = float(self._get(row, "rtol", 1e-10))
        atol = float(self._get(row, "atol", 1e10))
        fp_max_stepsize = float(self._get(row, "fp_max_stepsize", 5.0))
        max_stepsize_no_fp = float(self._get(row, "max_stepsize_no_fp", 100.0))

        bc_plasma_facing = self._get(row, "bc_plasma", "Robin - Surf. Rec. + Implantation")
        bc_rear = self._get(row, "bc_rear", "Neumann - no flux")

        calc_implant_str = self._get(row, "calc_implant", "Yes")
        calculate_implantation_params = str(calc_implant_str).lower().strip() != "no"

        # ── Sim ID (optional column, else 1-based row number) ──
        raw_sim_id = self._get(row, "sim_id", None)
        sim_id = int(raw_sim_id) if raw_sim_id is not None else (row_index + 1)

        # ── Atom view factor (optional, default 1.0) ──
        atom_view_factor = float(self._get(row, "atom_view_factor", 1.0))

        # ── Build config & Bin ──
        bin_config = BinConfiguration(
            rtol=rtol,
            atol=atol,
            fp_max_stepsize=fp_max_stepsize,
            max_stepsize_no_fp=max_stepsize_no_fp,
            bc_plasma_facing_surface=bc_plasma_facing,
            bc_rear_surface=bc_rear,
        )

        bin_obj = Bin(
            flux_id=flux_id,
            z_start=z_start,
            r_start=r_start,
            z_end=z_end,
            r_end=r_end,
            material=mat_obj,
            thickness=thickness,
            cu_thickness=cu_thickness,
            mode=mode,
            parent_bin_surf_area=parent_bin_surf_area,
            surface_area=surface_area,
            f_ion_flux_fraction=f_ion_flux_fraction,
            location=location,
            coolant_temp=coolant_temp,
            bin_configuration=bin_config,
            sim_id=sim_id,
            calculate_implantation_params=calculate_implantation_params,
            atom_view_factor=atom_view_factor,
        )
        return bin_obj
    
    def load_all_bins(self) -> BinCollection:
        """
        Load all bins from the CSV file.
        
        Returns:
            BinCollection containing all bins
        """
        bins = []
        for row_index, row in self.df.iterrows():
            bin_obj = self.load_bin_from_row(row, row_index)
            bins.append(bin_obj)
        print(f"✓ Successfully loaded {len(bins)} bins from CSV")
        return BinCollection(bins)
    
    def load_reactor(self) -> Reactor:
        """Load a complete reactor from the CSV file."""
        bin_collection = self.load_all_bins()
        return Reactor(bin_collection.bins)
    
    def get_summary(self) -> Dict[str, Any]:
        """Get summary statistics of the CSV data."""
        flux_col = self._csv_col("flux_id")
        mat_col  = self._csv_col("material")
        loc_col  = self._csv_col("location")
        mode_col = self._csv_col("mode")
        summary = {
            'total_rows': len(self.df),
            'unique_flux_ids': len(self.df[flux_col].unique()) if flux_col else 0,
            'materials': self.df[mat_col].value_counts().to_dict() if mat_col else {},
            'locations': self.df[loc_col].value_counts().to_dict() if loc_col else {},
            'modes': self.df[mode_col].value_counts().to_dict() if mode_col else {},
        }
        return summary
    
    def print_summary(self):
        """Print summary statistics of the CSV data."""
        summary = self.get_summary()
        
        print("\n=== CSV Data Summary ===")
        print(f"Total rows: {summary['total_rows']}")
        print(f"Unique flux IDs: {summary['unique_flux_ids']}")
        
        print("\nMaterials:")
        for material, count in summary['materials'].items():
            print(f"  {material}: {count}")
        
        print("\nLocations:")
        for location, count in summary['locations'].items():
            print(f"  {location}: {count}")
        
        print("\nModes:")
        for mode, count in summary['modes'].items():
            print(f"  {mode}: {count}")


def load_csv_reactor(csv_path: str, materials_csv_path: Optional[str] = None) -> Reactor:
    """
    Convenience function to load a reactor from CSV file.
    """
    if materials_csv_path is None:
        candidate = Path(csv_path).parent / "materials.csv"
        if candidate.exists():
            materials_csv_path = candidate
    loader = CSVBinLoader(csv_path, materials_csv_path=materials_csv_path)
    return loader.load_reactor()


# Example usage function
def example_usage():
    """Example of how to use the CSV bin loader."""
    csv_path = "input_files/input_table.csv"
    try:
        loader = CSVBinLoader(csv_path)
        loader.print_summary()
        reactor = loader.load_reactor()
        print(f"\n✓ Created {reactor}")
        print(f"  First wall bins: {len(reactor.first_wall_bins)}")
        print(f"  Divertor bins: {len(reactor.divertor_bins)}")
        if len(reactor.bins) > 0:
            example_bin = reactor.bins[0]
            print(f"\nExample bin: {example_bin}")
            print(f"  Configuration: rtol={example_bin.bin_configuration.rtol}")
            print(f"  Geometry: {example_bin.start_point} -> {example_bin.end_point}")
            print(f"  Length: {example_bin.length:.3f} m")
    except Exception as e:
        print(f"Error loading CSV reactor: {e}")


if __name__ == "__main__":
    example_usage()
