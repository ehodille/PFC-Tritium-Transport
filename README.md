# PFC-Tritium-Transport

The purpose of this repository is to obtain Hydrogen inventory estimations in the Plasma Facing Components in custom reactor and plasma scenarios. For now, it relies on HISP in order to run FESTIM simulations from which the inventory will be calculated. In the future this repository will be further developed to be able to work with different Hydrogen transport codes at user's preference. 

---

## Prerequisites

- **Operating System:** Linux (any modern distribution; tested on RHEL-based systems at ITER). macOS/Windows are not officially supported.
- **Conda distribution:** [Miniconda](https://docs.anaconda.com/miniconda/) is highly recommended. It provides a minimal conda installation and avoids the bloated dependency graph of full Anaconda, which can cause extremely slow or stalled environment solves.

---

## How to Install

### 1. Clone this repository
```bash
git clone --branch main https://github.com/iterorganization/PFC-Tritium-Transport.git
```

### 2. Create the conda environment
The `PFC-TT.yml` file included in this repository pins all required conda packages — including **FESTIM** and its heavy dependencies (FEniCSx, PETSc, MPI, etc.) — so that the solver resolves everything in a single pass.

```bash
conda env create -f PFC-TT.yml
conda activate PFC-TT
```

> **Important:** FESTIM must be solved together with all other conda packages at environment-creation time. Installing it separately afterwards (e.g. `conda install festim`) will typically hang indefinitely because the solver cannot reconcile FESTIM's dependency tree with an already-locked environment.

This step installs all core simulation dependencies, including:
- **FESTIM v2.0-beta.2**
- **FEniCS-DOLFINx v0.10.0** — finite element backend
- **FEniCS-Basix v0.10.0**, **FFCx v0.10.1**, **UFL v2025.2.1** — FEniCSx components
- **PETSc v3.24.5** + **petsc4py** — parallel linear algebra solvers
- **MPICH v4.3.2** + **mpi4py** — MPI support for parallel runs
- **HDF5 v1.14.6** — data I/O
- **pandas, numpy, scipy, matplotlib, jupyterlab** — data analysis and visualisation

### 3. Install HISP and h_transport_materials
We use a specific branch of HISP which will remain static until significant improvements are made.
To avoid FESTIM version conflicts, we install HISP **without dependencies**:
```bash
pip install --no-deps git+https://github.com/AdriaLlealS/hisp.git@main
pip install h_transport_materials
```

This will install:
- **HISP**
- **h_transport_materials**
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

#### Local (interactive)
Run `main.py` from the repository root. It will prompt for the input folder and the simulation ID:

```bash
python main.py
```

```
============================================================
  PFC-Tritium-Transport — local simulation runner
============================================================

Input folder path: my_input_folder
Simulation ID    : 1
```

- **Input folder path** — path to the folder containing `input_table.csv`, `scenario.py`, etc. Can be absolute or relative to the repository root.
- **Simulation ID** — integer ID of the bin to run, matching the `Sim. ID` column in `input_table.csv` (or 1-based row number if the column is absent).

Results are saved to `<input_folder>/results_<folder_name>/`.

#### Cluster (SLURM, ITER SCDCC)
- Use `run_on_cluster/slurm_folder_jobs.sh` to submit jobs to the cluster. This script is tailored for ITER's Science Division Computer Cluster (SCDCC) and includes site-specific module loads and paths.
- Example usages:

```bash
# Run all bins for a scenario (using input_files/ as default input folder)
sbatch run_on_cluster/slurm_folder_jobs.sh input_files

# Run specific bins for a scenario
sbatch run_on_cluster/slurm_folder_jobs.sh input_files "0-4, 10-12"
```

### Notes
- All data files (fluxes, heat loads, energies, etc.) are specified in `scenario.py`. You do not need to modify the runner or hardcode file paths elsewhere.
- The SLURM script is tailored for ITER's cluster; adapt it for other systems as needed.
- Column header names in `input_table.csv` are case-sensitive and must match exactly.
- Results are saved to a results folder inside your input folder (e.g. `my_input_folder/results_scenario/`).

---


