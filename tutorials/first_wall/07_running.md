# 07 — Running

**Previous:** [06 — Temperature Models](06_temperature_models.md) | **Next:** [08 — Sensitivity Scan](08_sensitivity_scan.md)

---

## Prerequisites

```bash
conda activate PFC-TT
export PFC_TT_PATH=$(pwd)          # run from repo root
```

## Running on the ITER cluster (SLURM)

The runner accepts an **absolute path**, a **relative path**, or a **bare folder name**
(searched inside `simulations/`). You do not have to be in the repo root:

```bash
# Option A — from repo root, using a relative path
bash run_on_cluster/slurm_folder_jobs.sh tutorials/first_wall/FW_example

# Option B — from anywhere, using an absolute path
bash /path/to/PFC-Tritium-Transport/run_on_cluster/slurm_folder_jobs.sh \
    /path/to/PFC-Tritium-Transport/tutorials/first_wall/FW_example
```

To submit a subset of sim IDs, append the range after the folder path:

```bash
bash run_on_cluster/slurm_folder_jobs.sh tutorials/first_wall/FW_example "1-2"
```

The script:
1. Reads `input_table.csv` to find all sim IDs.
2. Submits one `sbatch` job per sim ID (`--time=300h`, `--mem=1gb`, `--ntasks=1`).
3. Writes logs to `<input_folder>/logs/sim_<id>_<jobid>.out/.err`.

## Running locally (single bin, no SLURM)

The `run_new_csv_bin.py` script also accepts absolute paths:

```bash
# From repo root
python run_on_cluster/run_new_csv_bin.py \
    1 \                                           # Sim. ID
    tutorials/first_wall/FW_example \             # input folder (relative or absolute)
    scenario_reference \                           # scenario module name (no .py)
    tutorials/first_wall/FW_example/input_table.csv \
    --input-dir tutorials/first_wall/FW_example

# From anywhere (absolute paths)
FW=/path/to/PFC-Tritium-Transport/tutorials/first_wall/FW_example
python /path/to/PFC-Tritium-Transport/run_on_cluster/run_new_csv_bin.py \
    1 "$FW" scenario_reference "$FW/input_table.csv" --input-dir "$FW"
```

> Running all 3 bins locally will take several hours each for the full
> 35-pulse scenario. For a quick test, reduce `nb_pulses` in `scenario_reference.py`.

## Monitoring progress

```bash
python check_progress.py tutorials/first_wall/FW_example
```

Output:
```
=== Progress: tutorials/first_wall/FW_example ===
Completed:  1 / 3
Running:    1 / 3
Failed:     0 / 3

Sim ID | Elapsed | Remaining | Status
     1 |   2h14m |    ~8h32m | running
     2 |       — |         — | queued
     3 | finished|         — | completed
```

See [wiki/Running-on-Cluster.md](../../wiki/Running-on-Cluster.md) for full monitoring details.

## Output folder structure

Results are written to `tutorials/first_wall/FW_example/results_FW_example/`:

```
results_FW_example/
├── sim_1_flux_0_w_wetted.json
├── sim_2_flux_1_w_wetted.json
└── sim_3_flux_2_w_shadowed.json

profiles_FW_example/
├── sim_1_flux_0_w_wetted_profiles.json
└── ...

checkpoints/
├── checkpoint_sim_1_flux_0.bp
└── ...
```

## Expected results (indicative)

After completing, examine the JSON output for each bin:

| Sim. ID | Bin | mode | Peak T-inventory (D+T at./m²) | Notes |
|---------|-----|------|-------------------------------|-------|
| 1 | 0 | wetted | ~10²² | moderate retention |
| 2 | 1 | wetted | ~2×10²² | highest flux → highest inventory |
| 3 | 2 | shadowed | ~10¹⁸ | negligible flux, near-zero retention |

> These are order-of-magnitude estimates for context. Actual values depend on
> the full scenario timing and tritium fraction.

**Physically meaningful signs:**
- Inventory grows during pulses and partially releases during baking.
- Permeation flux (at rear surface, if Neumann BC is used) is much lower than
  implanted flux — most tritium is trapped.
- The shadowed bin (ID 3) shows nearly flat retention over time.

## Common issues

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `plasma_data_handling is None` | Wrong path to `.dat` files | Check `_data_dir` in scenario; run `python -c "from pathlib import Path; print(Path('tutorials/first_wall/FW_example/scenario_reference.py').resolve().parent.parent / 'plasma_data')"` |
| `KeyError: 'FP_ALT'` | `PlasmaDataHandling` missing the key | Ensure `pulse_type_to_data` has `"FP_ALT"` matching the `Pulse(pulse_type=...)` |
| Solver divergence in `.err` log | Time step too large during ramp | Reduce `FP max. stepsize` in `input_table.csv` |
| No `.json` output but no error | Job still running | Check `squeue` or `check_progress.py` |

---

**Next:** [08 — Sensitivity Scan](08_sensitivity_scan.md)
