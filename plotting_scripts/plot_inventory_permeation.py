"""
Plot tritium inventory and permeation flux over time for a single simulation bin.

Subplots:
  1. Total retention, trap inventory, and mobile inventory vs time [atoms/m²]
  2. Plasma-facing surface temperature vs time [K]
  3. Instantaneous permeation flux (left y-axis) and cumulative permeation
     (right y-axis) vs time.  Uses T_outlet_flux for permeation.

Data source:
  <RESULTS_FOLDER>/sim_<SIM_ID>_*.json
  Keys used:
    "t"                 – list[float], time in seconds
    "T"                 – {"data": [...]}, mobile T inventory (atoms/m²)
    "trap1_T", "trap2_T", ...
                        – {"data": [...]}, trapped T inventory (atoms/m²)
    "T_outlet_flux"     – {"data": [...]}, rear surface T flux (atoms/m²/s)
    "temperature_at_x0" – list[float], plasma-facing surface temperature (K)
  Metadata keys: "sim_id", "mode", "material", "location"

Output is saved to <scenario_dir>/plots/  (sibling of the results folder).
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ---------------------------------------------------------------------------
# CONFIG – edit these two variables before running
# ---------------------------------------------------------------------------
RESULTS_FOLDER = Path(
    "../tutorials/first_wall/FW_example/results_FW_example"
)  # path to the results folder (absolute or relative to this script)
SIM_ID = 1  # integer – the first number in the result filename, e.g. 1 for sim_1_*
# ---------------------------------------------------------------------------

FIGSIZE = (12, 9)
TIME_UNIT = "days"  # display unit for time axis; seconds are converted below


def resolve_folder(p: Path) -> Path:
    """Resolve path relative to this script file if not absolute."""
    return p if p.is_absolute() else (Path(__file__).parent / p).resolve()


def seconds_to_display(t_s: np.ndarray) -> np.ndarray:
    if TIME_UNIT == "days":
        return t_s / 86400.0
    if TIME_UNIT == "hours":
        return t_s / 3600.0
    return t_s  # seconds


def extract_array(obj) -> np.ndarray:
    """Return a float array from either {"data": [...]} or a plain list."""
    if isinstance(obj, dict):
        return np.array(obj.get("data", []), dtype=float)
    return np.array(obj, dtype=float)


def find_json(results_folder: Path, sim_id: int) -> Path:
    """Glob for sim_<sim_id>_*.json and return the first match."""
    matches = sorted(results_folder.glob(f"sim_{sim_id}_*.json"))
    if not matches:
        raise FileNotFoundError(
            f"No file matching 'sim_{sim_id}_*.json' in {results_folder.resolve()}"
        )
    if len(matches) > 1:
        print(f"[WARN] Multiple matches for sim_id={sim_id}: {[m.name for m in matches]}. Using {matches[0].name}")
    return matches[0]


def main():
    results_folder = resolve_folder(RESULTS_FOLDER)
    json_path = find_json(results_folder, SIM_ID)
    print(f"Loading: {json_path.name}")

    with open(json_path) as f:
        data = json.load(f)

    if "t" not in data:
        raise KeyError("JSON is missing the required 't' (time) array.")

    # ---- Time ---------------------------------------------------------------
    t_s = np.array(data["t"], dtype=float)
    t = seconds_to_display(t_s)

    # ---- Inventory arrays ---------------------------------------------------
    mobile_T = extract_array(data.get("T", []))

    trap_keys = sorted(k for k in data if k.startswith("trap") and k.endswith("_T"))
    traps = {k: extract_array(data[k]) for k in trap_keys}
    trap_T_total = sum(traps.values()) if traps else np.zeros_like(mobile_T)

    total_retention = mobile_T + trap_T_total

    # ---- Surface temperature ------------------------------------------------
    temp_at_x0 = extract_array(data.get("temperature_at_x0", []))

    # ---- Permeation flux (rear/outlet) -------------------------------------
    # T_outlet_flux: rear surface flux, atoms/m²/s  (positive = flux away from plasma)
    # TODO: verify sign convention – some codes may export negative outlet flux
    T_outlet_flux = extract_array(data.get("T_outlet_flux", []))
    has_permeation = T_outlet_flux.size > 0

    if has_permeation:
        # Cumulative permeation via numerical integration
        cumulative_T = np.cumsum(T_outlet_flux * np.gradient(t_s))
        # cumulative_T units: atoms/m² (flux [atoms/m²/s] × dt [s])

    # ---- Metadata for title -------------------------------------------------
    sim_id_val = data.get("sim_id", SIM_ID)
    mode       = data.get("mode", "unknown")
    material   = data.get("material", "?")
    location   = data.get("location", "")
    title = f"{location} sim {sim_id_val}, {mode}, {material}  [{json_path.stem}]"

    # ---- Plot ---------------------------------------------------------------
    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, sharex=True, figsize=FIGSIZE,
        gridspec_kw={"height_ratios": [3, 1, 2]},
    )

    # Subplot 1: inventory
    color_mobile = "tab:gray"
    ax1.plot(t, total_retention, color="black", linewidth=2, label="Total T inventory")
    for trap_key, trap_arr in traps.items():
        ax1.plot(t, trap_arr, linewidth=1.2, linestyle="-", label=trap_key)

    ax1.set_ylabel("T Inventory (atoms/m²)")
    ax1.set_title(title)

    ax1_mobile = ax1.twinx()
    ax1_mobile.plot(t, mobile_T, color=color_mobile, linewidth=0.8,
                    linestyle="--", label="Mobile T")
    ax1_mobile.set_ylabel("Mobile T (atoms/m²)", color=color_mobile)
    ax1_mobile.tick_params(axis="y", labelcolor=color_mobile)

    # Combined legend
    lines_l, labels_l = ax1.get_legend_handles_labels()
    lines_r, labels_r = ax1_mobile.get_legend_handles_labels()
    ax1.legend(lines_l + lines_r, labels_l + labels_r, loc="upper left", fontsize=9)

    ax1.grid(which="both", linestyle="--", linewidth=0.5)

    # Subplot 2: surface temperature
    if temp_at_x0.size > 0:
        ax2.plot(t, temp_at_x0, color="tab:red", linewidth=1.2,
                 label="Surface temperature")
        ax2.set_ylabel("T$_{surface}$ (K)")
        ax2.legend(loc="upper left")
    else:
        ax2.text(0.5, 0.5, "No temperature_at_x0 data",
                 transform=ax2.transAxes, ha="center", va="center", color="gray")
        ax2.set_ylabel("T$_{surface}$ (K)")
    ax2.grid(which="both", linestyle="--", linewidth=0.5)

    # Subplot 3: permeation (twin axes)
    if has_permeation:
        color_inst = "tab:purple"
        color_cum  = "tab:green"

        ax3.plot(t, T_outlet_flux, color=color_inst, linewidth=1.2,
                 label="Inst. permeation flux")
        ax3.set_ylabel("Permeation flux\n(atoms/m²/s)", color=color_inst)
        ax3.tick_params(axis="y", labelcolor=color_inst)

        ax3_right = ax3.twinx()
        ax3_right.plot(t, cumulative_T, color=color_cum, linewidth=1.5,
                       linestyle="--", label="Cumul. permeation")
        ax3_right.set_ylabel("Cumulative permeation\n(atoms/m²)", color=color_cum)
        ax3_right.tick_params(axis="y", labelcolor=color_cum)

        # Combined legend
        lines_l, labels_l = ax3.get_legend_handles_labels()
        lines_r, labels_r = ax3_right.get_legend_handles_labels()
        ax3.legend(lines_l + lines_r, labels_l + labels_r,
                   loc="upper left", fontsize=9)
    else:
        ax3.text(0.5, 0.5, "No T_outlet_flux data available\n(permeation not exported)",
                 transform=ax3.transAxes, ha="center", va="center",
                 fontsize=10, color="gray")
        ax3.set_ylabel("Permeation flux (atoms/m²/s)")

    ax3.set_xlabel(f"Time ({TIME_UNIT})")
    ax3.grid(which="both", linestyle="--", linewidth=0.5)

    plt.tight_layout()

    # Save to <scenario_dir>/plots/
    plots_dir = results_folder.parent / "plots" / f"sim_{SIM_ID}"
    plots_dir.mkdir(parents=True, exist_ok=True)
    out_path = plots_dir / f"{json_path.stem}_inventory_permeation.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.show()


if __name__ == "__main__":
    main()
