# PFC-Tritium-Transport

The purpose of this repository is to obtain Hydrogen inventory estimations in the Plasma Facing Components in custom reactor and plasma scenarios. For now, it relies on HISP in order to run FESTIM simulations from which the inventory will be calculated. In the future this repository will be further developed to be able to work with different Hydrogen transport codes at user's preference. 

---

## How to Install

### 1. Clone this repository
```bash
git clone --branch main https://github.com/iterorganization/PFC-Tritium-Transport.git
```

### 2. Recreate the Python environment
We use **conda** for reproducibility. Make sure the file `PFC-TT.yml` was correctly downloaded from this repository and run:
```bash
conda env create -f PFC-TT.yml
conda activate PFC-TT
```

---

### 3. Install FESTIM2 and custom HISP version
We use a specific branch of HISP which will remain static until significant improvements are made. 
To avoid FESTIM version conflicts, we install HISP **without dependencies**:
```bash
conda install -c conda-forge 'festim=2.0b2.post2'
pip install --no-deps git+https://github.com/AdriaLlealS/hisp.git@main
pip install h_transport_materials
```
This will install:
- **FESTIM v2.0-beta.2**
- **dolfinx v0.10.0**
- **HISP**
- All required dependencies

---

## How to Run (Input Folder Methodology)

All simulations are now launched from a dedicated input folder containing the following files:

- **input_table.csv** (required):
  - Describes all bins to simulate. Required columns (case-sensitive):
    - `Bin number`, `Z_start (m)`, `R_start (m)`, `Z_end (m)`, `R_end (m)`,
    - `Material`, `Thickness (m)`, `Cu thickness (m)`, `mode`,
    - `S. Area parent bin (m^2)`, `Surface area (m^2)`, `f (ion flux scaling factor)`, `location`.
    - Optional: `BC Plasma Facing Surface`, `BC rear surface`, `rtol`, `atol`, `FP max. stepsize (s)`, `Max. stepsize no FP (s)`, etc.
- **materials.py** (required):
  - Python file defining all materials referenced in the input table. Each material should be a Python object or dictionary.
- **scenario.py** (required):
  - Python file defining the scenario (pulse sequence, campaign structure, etc). **All data file paths for fluxes, heat loads, particle energies, etc. are specified here.** This is the only place you need to set these paths.
- **mesh.py** (optional but recommended):
  - Python file describing the mesh for each bin or group. If omitted, a default mesh is generated from the input table.
- **temperature_models.py** (optional):
  - Python file with custom temperature models (functions) for one or more materials/bins. If omitted, default temperature handling is used.

### Workflow Summary

1. **Prepare your input folder** with the files above. Only `input_table.csv`, `materials.py`, and `scenario.py` are strictly required.
2. **Edit `scenario.py`** to reference the correct data files for fluxes, heat loads, particle energies, etc. (e.g. `Binned_Flux_Data.dat`, `ICWC_data.dat`, etc.).
3. **(Optional) Add `mesh.py` and/or `temperature_models.py`** for custom mesh or temperature logic.
4. **Run the simulation** using one of the runners below, pointing to your input folder.

### Runners: Cluster and Local

#### Cluster (SLURM, ITER SCDCC)
- Use `run_on_cluster/slurm_folder_jobs.sh` to submit jobs to the cluster. This script is tailored for ITER's Science Division Computer Cluster (SCDCC) and includes site-specific module loads and paths.
- Example usages:

```bash
# Run all bins for a scenario (using input_files/ as default input folder)
sbatch run_on_cluster/slurm_folder_jobs.sh input_files

# Run specific bins for a scenario
sbatch run_on_cluster/slurm_folder_jobs.sh input_files "0-4, 10-12"
```

#### Local (Single Bin, Any System)
- Use `run_on_cluster/run_bin_from_folder.py` to run a single bin locally from an input folder. This is useful for testing or debugging without SLURM.
- Example usage:

```bash
python run_on_cluster/run_bin_from_folder.py my_input_folder 1
```
- This will run the simulation for bin ID 1 (first row in `input_table.csv`) using the files in `my_input_folder/`.

### Notes
- All data files (fluxes, heat loads, energies, etc.) are specified in `scenario.py`. You do not need to modify the runner or hardcode file paths elsewhere.
- The SLURM script is tailored for ITER's cluster; adapt it for other systems as needed.
- Column header names in `input_table.csv` are case-sensitive and must match exactly.
- Results are saved to a results folder inside your input folder (e.g. `my_input_folder/results_scenario/`).

---


