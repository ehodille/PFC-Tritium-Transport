# 08 — Sensitivity Scan

**Previous:** [07 — Running](07_running.md) | **Top:** [README](README.md)

---

Two common scan types are demonstrated below:
1. **Scenario scan** — swap `scenario_reference.py` for `scenario_highflux.py`
2. **Material scan** — change `W` to `W_hightraps` in `input_table.csv`

## 1. Scenario scan (reference vs. high-flux)

### Setup

```bash
# Duplicate the input folder
cp -r tutorials/first_wall/FW_example tutorials/first_wall/FW_example_highflux

# Swap the scenario file
rm tutorials/first_wall/FW_example_highflux/scenario_reference.py
cp tutorials/first_wall/scenarios/scenario_highflux.py \
   tutorials/first_wall/FW_example_highflux/
```

### Run

```bash
# Reference (already done, or re-run):
bash run_on_cluster/slurm_folder_jobs.sh tutorials/first_wall/FW_example

# High-flux variant:
bash run_on_cluster/slurm_folder_jobs.sh tutorials/first_wall/FW_example_highflux
```

Results land in:
- `tutorials/first_wall/FW_example/results_FW_example/`
- `tutorials/first_wall/FW_example_highflux/results_FW_example_highflux/`

### Compare

```python
import json, pathlib

def load_inventory(results_dir, sim_id, flux_id, material, mode):
    fname = f"sim_{sim_id}_flux_{flux_id}_{material}_{mode}.json"
    data = json.loads((pathlib.Path(results_dir) / fname).read_text())
    # Total inventory = sum of all trapped species at final time step
    return sum(v["data"][-1] for k, v in data.items()
               if k.startswith("trap") and isinstance(v, dict))

ref_inv  = load_inventory("...results_FW_example",          1, 0, "w", "wetted")
hf_inv   = load_inventory("...results_FW_example_highflux", 1, 0, "w", "wetted")
print(f"Ratio high-flux / reference: {hf_inv / ref_inv:.2f}")
```

The high-flux scenario applies `flux_scaling=2.0` and `heat_scaling=1.5` during
the 20 full-power pulses, so you should see roughly 1.5–2× higher retention
(sub-linear because trapping sites saturate).

---

## 2. Material parameter scan (trap density)

### Setup

In `tutorials/first_wall/FW_example/input_table.csv`, change the `Material`
column for some bins:

```
# Reference: all bins use "W"
1,0,...,W,...
2,1,...,W,...

# Sensitivity: replace with "W_hightraps" (5× trap density)
1,0,...,W_hightraps,...
2,1,...,W_hightraps,...
```

`W_hightraps` is already defined in `materials.csv` — no other file changes needed.

Run as before:

```bash
bash run_on_cluster/slurm_folder_jobs.sh tutorials/first_wall/FW_example
```

### Why define multiple material variants in materials.csv?

This pattern is the cleanest way to scan material parameters:

```
Material_name, W, , Material_name, W_hightraps, , Material_name, W_frauenfelder
Mat_density,   6.3382E+28, , Mat_density, 6.3382E+28, , Mat_density, 6.3382E+28
D0,            2.06E-07,   , D0,          2.06E-07,   , D0,          4.1E-07
...
```

Then in `input_table.csv` you can have rows for the **same bin** with different
material names:

```
Sim. ID, Bin number, ..., Material, ...
1,       0,          ..., W,        ...
4,       0,          ..., W_hightraps, ...
5,       0,          ..., W_frauenfelder, ...
```

Each row gets a unique `Sim. ID` → separate output JSON → easy comparison.

---

## General scan workflow

```
1. Choose the parameter to vary.
2. Either:
   a. Duplicate the input folder and change the scenario file (scenario scan), or
   b. Add material variants to materials.csv and add rows to input_table.csv (material scan).
3. Run with slurm_folder_jobs.sh.
4. Load the resulting JSON files and compare the quantity of interest.
```

> Results from different runs are independent and can be compared
> post-hoc without re-running; results JSONs are self-contained (they include
> all metadata needed to identify the case).

---

**Top:** [README](README.md) | **Wiki:** [wiki/Home.md](../../wiki/Home.md)
