# 09 — Plotting Results

**Previous:** [08 — Sensitivity Scan](08_sensitivity_scan.md) | **Top:** [README](README.md)

---

Three standalone plotting scripts live in `plotting_scripts/` at the root of the
repository.  Each script has two (or three) configuration variables at the top
that you edit before running — no command-line arguments needed.

| Script | What it plots |
|--------|---------------|
| `plot_inventory.py` | Tritium inventory over time (retention, traps, mobile) + surface temperature |
| `plot_inventory_permeation.py` | Same as above, plus a permeation flux panel |
| `plot_profile.py` | Concentration profile through the material at one timestep |

Plots are saved automatically to `<case_folder>/plots/sim_<N>/`.

---

## Configuration variables

Every script has a `CONFIG` block near the top:

```python
RESULTS_FOLDER = Path("../tutorials/first_wall/FW_example/results_FW_example")
# path to the results folder (absolute or relative to the script file)

SIM_ID = 1
# integer – the first number in the result filename, e.g. 1 for sim_1_*
```

`plot_profile.py` has one extra variable:

```python
TIMESTEP_INDEX = -1
# integer index into the stored timesteps (0 = first, -1 = last)
```

---

## How to run

Run any script directly with Python from any working directory:

```bash
python plotting_scripts/plot_inventory.py
python plotting_scripts/plot_inventory_permeation.py
python plotting_scripts/plot_profile.py
```

To generate plots for all simulations in a case at once:

```python
import importlib.util

scripts = [
    "plotting_scripts/plot_inventory.py",
    "plotting_scripts/plot_inventory_permeation.py",
    "plotting_scripts/plot_profile.py",
]

for sim_id in [1, 2, 3]:
    for script_path in scripts:
        spec = importlib.util.spec_from_file_location("mod", script_path)
        mod  = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        mod.SIM_ID = sim_id
        mod.main()
```

---

## Script details

### `plot_inventory.py` — inventory over time

Two subplots sharing the time axis:

- **Top** — tritium inventory (atoms/m²):
  - Total retention (black, thick solid)
  - One line per trap (`trap1_T`, `trap2_T`, …) on the left y-axis
  - Mobile T concentration on an independent right y-axis (gray, thin dashed)
- **Bottom** — plasma-facing surface temperature (K)

<img width="1781" height="880" alt="sim_1_flux_0_w_wetted_inventory" src="https://github.com/user-attachments/assets/51d43219-27d9-46cc-8f7b-9b62d3306b5c" />

---

### `plot_inventory_permeation.py` — inventory + permeation

Same as above but with a third subplot:

- **Top** — tritium inventory (same as above)
- **Middle** — surface temperature (K)
- **Bottom** — permeation through the rear surface:
  - Instantaneous flux in atoms/m²/s (left y-axis, purple)
  - Cumulative permeation in atoms/m² (right y-axis, green dashed)

<img width="1781" height="1331" alt="sim_1_flux_0_w_wetted_inventory_permeation" src="https://github.com/user-attachments/assets/7ea63a25-dda7-4931-bbc5-86c3e19f66ac" />

---

### `plot_profile.py` — concentration profile at one timestep

Single figure with twin y-axes:

- **x-axis** — depth into the material (mm)
- **Left y-axis** — mobile T concentration (m⁻³, blue)
- **Right y-axis** — total trapped T concentration (m⁻³, red)

The figure title shows the bin label and the physical time corresponding to
`TIMESTEP_INDEX`.

> **Note:** profile data is read from `profiles_<case>/sim_<N>_*_profiles.json`,
> which is derived automatically from `RESULTS_FOLDER` by replacing the
> `results_` prefix with `profiles_`.

<img width="1483" height="882" alt="sim_1_flux_0_w_wetted_profiles_t244" src="https://github.com/user-attachments/assets/1dde155d-34dd-4cbf-a10f-07740d6ea002" />


---

**Previous:** [08 — Sensitivity Scan](08_sensitivity_scan.md) | **Top:** [README](README.md) | **Wiki:** [Results reference](../../wiki/Results.md)
