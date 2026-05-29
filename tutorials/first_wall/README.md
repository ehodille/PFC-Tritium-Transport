# Tutorial: First Wall Tritium Transport

A self-contained worked example simulating tritium retention in three ITER first-wall bins over a 35-pulse scenario with a baking event.

**What you will learn:**
- How the input folder structure works end to end
- How to define a multi-phase scenario with two plasma data files
- How to run the simulation and monitor progress
- How to perform a simple parameter scan

## Prerequisites

Environment and installation: see [wiki/Home.md](../../wiki/Home.md).

```bash
conda activate PFC-TT
export PFC_TT_PATH=$(pwd)   # must be run from the repo root
```

## Folder structure

```
tutorials/first_wall/
├── plasma_data/
│   ├── fp_pulse.dat                 ← standard FP plasma data (3 FW bins)
│   └── alternative_fp_pulse.dat    ← alternative magnetic config data
├── FW_example/                      ← the folder you pass to the runner
│   ├── input_table.csv
│   ├── materials.csv
│   ├── mesh.py
│   ├── temperature_models.py
│   └── scenario_reference.py       ← active scenario (copy/replace to change)
└── scenarios/
    └── scenario_highflux.py        ← sensitivity variant (not run by default)
```

## Walkthrough (read in order)

| Step | File | Topic |
|------|------|-------|
| 1 | [01_plasma_data.md](01_plasma_data.md) | Plasma flux data files |
| 2 | [02_input_table.md](02_input_table.md) | Bin definitions (`input_table.csv`) |
| 3 | [03_materials.md](03_materials.md) | Material parameters (`materials.csv`) |
| 4 | [04_scenario.md](04_scenario.md) | Scenario and pulse sequence |
| 5 | [05_mesh.md](05_mesh.md) | Spatial mesh (`mesh.py`) |
| 6 | [06_temperature_models.md](06_temperature_models.md) | Custom temperature model |
| 7 | [07_running.md](07_running.md) | Running and monitoring |
| 8 | [08_sensitivity_scan.md](08_sensitivity_scan.md) | Parameter scan |
| 9 | [09_plotting_results.md](09_plotting_results.md) | Plotting results |

## Quick start (if you already understand the code)

```bash
# From repo root
bash run_on_cluster/slurm_folder_jobs.sh tutorials/first_wall/FW_example
```

Results appear in `tutorials/first_wall/FW_example/results_FW_example/`.
