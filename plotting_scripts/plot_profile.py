"""
Plot a concentration profile through the material thickness for one bin at one
timestep.

Profile data source:
  <profiles_folder>/sim_<SIM_ID>_*_profiles.json
  The profiles folder is derived from RESULTS_FOLDER by replacing the
  'results_' prefix with 'profiles_' in the folder name.

  Keys used:
    "T_profile"            – {"x": [...], "t": [...], "data": [[...], ...]},
                             mobile T concentration (m⁻³);
                             x in metres, t in seconds
    "trap1_T_profile",
    "trap2_T_profile", ... – same structure, but x array is duplicated
                             (every-other point is taken; see note below)
    "D_profile"            – present but NOT plotted here; only T is shown
                             # TODO: verify – add D if needed

  Note on trap x-arrays: the existing codebase (plot_profiles.py) shows that
  trap profile arrays are stored with twice as many x-points as the mobile
  profile.  A stride-2 slice (arr[::2]) is used to de-duplicate.
  If this does not produce the right spatial grid, adjust the dedup logic.
  # TODO: verify deduplication for different material types

x-axis: depth into material (mm, converted from metres)
y-axis: concentration (m⁻³)

Mobile T is plotted on the left axis; total trapped T on the right axis.
The figure title shows the bin name and the physical time corresponding to
TIMESTEP_INDEX.

Output is saved to <scenario_dir>/plots/  (sibling of the results folder).
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import re

# ---------------------------------------------------------------------------
# CONFIG – edit these three variables before running
# ---------------------------------------------------------------------------
RESULTS_FOLDER = Path(
    "../tutorials/first_wall/FW_example/results_FW_example"
)  # path to the results folder (absolute or relative to this script)
SIM_ID = 1  # integer – the first number in the result filename, e.g. 1 for sim_1_*
TIMESTEP_INDEX = -1  # integer index into the stored timesteps (0 = first, -1 = last)
# ---------------------------------------------------------------------------


def resolve_folder(p: Path) -> Path:
    """Resolve path relative to this script file if not absolute."""
    return p if p.is_absolute() else (Path(__file__).parent / p).resolve()


def find_json(folder: Path, sim_id: int, suffix: str = "") -> Path:
    """Glob for sim_<sim_id>_*{suffix}.json and return the first match."""
    pattern = f"sim_{sim_id}_*{suffix}.json"
    matches = sorted(folder.glob(pattern))
    if not matches:
        raise FileNotFoundError(
            f"No file matching '{pattern}' in {folder.resolve()}"
        )
    if len(matches) > 1:
        print(f"[WARN] Multiple matches: {[m.name for m in matches]}. Using {matches[0].name}")
    return matches[0]


def deduplicate_x(x_arr: np.ndarray, ref_len: int) -> tuple[np.ndarray, np.ndarray]:
    """Return a de-duplicated x array and an index mask for the data.

    Trap profile x arrays are stored with ~2× spatial points.  Taking every
    other value matches the mobile grid.
    # TODO: verify – if the plot looks wrong, inspect the x arrays directly.
    """
    if len(x_arr) == ref_len:
        return x_arr, np.arange(ref_len)
    # Stride-2 slice
    x_dedup = x_arr[::2]
    idx     = np.arange(0, len(x_arr), 2)
    # Pad or trim to match ref_len
    if len(x_dedup) < ref_len:
        x_dedup = np.pad(x_dedup, (0, ref_len - len(x_dedup)), mode="edge")
        idx     = np.pad(idx,     (0, ref_len - len(idx)),     mode="edge")
    elif len(x_dedup) > ref_len:
        x_dedup = x_dedup[:ref_len]
        idx     = idx[:ref_len]
    return x_dedup, idx


def main():
    # Resolve RESULTS_FOLDER relative to the script file if it is not absolute
    results_folder  = resolve_folder(RESULTS_FOLDER)
    profiles_folder = results_folder.parent / re.sub(r"^results", "profiles", results_folder.name)
    profiles_path   = find_json(profiles_folder, SIM_ID, suffix="_profiles")

    if not profiles_path.exists():
        raise FileNotFoundError(
            f"Profiles JSON not found: {profiles_path.resolve()}\n"
            "Check that RESULTS_FOLDER and SIM_ID are correct, and that the "
            "profiles folder exists alongside the results folder."
        )

    with open(profiles_path) as f:
        pdata = json.load(f)

    # ---- Mobile T profile ---------------------------------------------------
    if "T_profile" not in pdata:
        raise KeyError(
            f"'T_profile' key not found in {profiles_path.name}.\n"
            f"Available keys: {list(pdata.keys())}"
        )

    T_prof = pdata["T_profile"]
    x_m    = np.array(T_prof["x"], dtype=float)   # metres
    t_s    = np.array(T_prof["t"], dtype=float)    # seconds
    n_t    = len(t_s)
    n_x    = len(x_m)

    # Resolve negative index
    idx = TIMESTEP_INDEX if TIMESTEP_INDEX >= 0 else n_t + TIMESTEP_INDEX
    if not (0 <= idx < n_t):
        raise IndexError(
            f"TIMESTEP_INDEX={TIMESTEP_INDEX} resolves to {idx}, "
            f"but there are only {n_t} timesteps (0 … {n_t - 1})."
        )

    mobile_conc = np.array(T_prof["data"][idx], dtype=float)
    if len(mobile_conc) != n_x:
        raise ValueError(
            f"Spatial dimension mismatch: x has {n_x} points but data row has "
            f"{len(mobile_conc)} points.  # TODO: verify data reshaping needed."
        )

    t_val_s    = t_s[idx]
    t_val_days = t_val_s / 86400.0
    x_mm       = x_m * 1e3   # convert to mm for display

    # ---- Trap T profiles ----------------------------------------------------
    trap_prof_keys = sorted(
        k for k in pdata if re.match(r"trap\d+_T_profile", k)
    )

    total_trapped = np.zeros(n_x, dtype=float)
    trap_found    = False

    for trap_key in trap_prof_keys:
        td = pdata[trap_key]
        tx = np.array(td["x"], dtype=float)
        ta = np.array(td["data"][idx], dtype=float)

        if len(ta) != n_x:
            # De-duplicate x dimension (trap arrays are stored with 2× x points)
            _, dedup_idx = deduplicate_x(tx, n_x)
            if len(ta) >= len(dedup_idx):
                ta = ta[dedup_idx]
            else:
                print(f"[WARN] Skipping {trap_key}: cannot reconcile spatial dimension.")
                continue

        total_trapped += ta
        trap_found = True

    # ---- Metadata from the results JSON (optional, for title) ---------------
    results_json = find_json(results_folder, SIM_ID)
    sim_label = profiles_path.stem  # fallback
    with open(results_json) as f:
        rdata = json.load(f)
    sim_id_val = rdata.get("sim_id", SIM_ID)
    mode       = rdata.get("mode", "?")
    material   = rdata.get("material", "?")
    location   = rdata.get("location", "")
    sim_label = f"{location} sim {sim_id_val}, {mode}, {material}"

    title = (
        f"{sim_label}"
        f"  [t = {t_val_s:.2f} s  /  {t_val_days:.3f} days,"
        f" timestep {idx}/{n_t - 1}]"
    )

    # ---- Plot ---------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(10, 6))

    color_mobile = "tab:blue"
    color_trap   = "tab:red"

    ax1.plot(x_mm, mobile_conc, color=color_mobile, linewidth=2,
             label="Mobile T")
    ax1.set_xlabel("Depth (mm)")
    ax1.set_ylabel("Mobile T concentration (m⁻³)", color=color_mobile)
    ax1.tick_params(axis="y", labelcolor=color_mobile)
    ax1.grid(which="both", linestyle="--", linewidth=0.5, alpha=0.6)

    if trap_found:
        ax2 = ax1.twinx()
        ax2.plot(x_mm, total_trapped, color=color_trap, linewidth=2,
                 label="Total trapped T")
        ax2.set_ylabel("Trapped T concentration (m⁻³)", color=color_trap)
        ax2.tick_params(axis="y", labelcolor=color_trap)

        # Combined legend
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc="best", fontsize=9)
    else:
        ax1.legend(loc="best")
        print("[WARN] No trap profile data found; only mobile T is plotted.")

    ax1.set_title(title, fontsize=10)
    plt.tight_layout()

    # Save to <scenario_dir>/plots/
    plots_dir = results_folder.parent / "plots" / f"sim_{SIM_ID}"
    plots_dir.mkdir(parents=True, exist_ok=True)
    out_path = plots_dir / f"{profiles_path.stem}_t{idx}.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.show()


if __name__ == "__main__":
    main()
