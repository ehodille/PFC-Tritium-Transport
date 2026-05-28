"""
scenario_highflux.py — Tutorial: First Wall (high-flux sensitivity variant)
============================================================================
Identical to scenario_reference.py except that the full-power phases use
2× particle flux and 1.5× heat load — representing a higher-performance
plasma or a more optimistic flux extrapolation.

To use this scenario:
    cp tutorials/first_wall/scenarios/scenario_highflux.py \\
       tutorials/first_wall/FW_example_highflux/
    rm tutorials/first_wall/FW_example_highflux/scenario_reference.py
    bash run_on_cluster/slurm_folder_jobs.sh tutorials/first_wall/FW_example_highflux

See 08_sensitivity_scan.md for the full comparison workflow.
"""

from pathlib import Path
from scenario import Scenario, Pulse
import pandas as pd
from plasma_data_handling import PlasmaDataHandling

# ── Plasma data ───────────────────────────────────────────────────────────────
_data_dir = Path(__file__).resolve().parent.parent / "plasma_data"

plasma_data_handling = PlasmaDataHandling(
    pulse_type_to_data={
        "FP":     pd.read_csv(_data_dir / "fp_pulse.dat"),
        "FP_ALT": pd.read_csv(_data_dir / "alternative_fp_pulse.dat"),
    },
)

# ── Pulses ────────────────────────────────────────────────────────────────────
# Phase 1: same as reference (commissioning phase — low power)
phase1_fp = Pulse(
    pulse_type="FP",
    nb_pulses=10,
    ramp_up=400,
    steady_state=600,
    ramp_down=400,
    waiting=3600,
    tritium_fraction=0.5,
    heat_scaling=0.5,
    flux_scaling=0.5,
)

# Phase 2: HIGH-FLUX — 2× particle flux, 1.5× heat load
phase2_fp = Pulse(
    pulse_type="FP",
    nb_pulses=20,
    ramp_up=400,
    steady_state=600,
    ramp_down=400,
    waiting=3600,
    tritium_fraction=0.5,
    heat_scaling=1.5,   # ← changed from 1.0
    flux_scaling=2.0,   # ← changed from 1.0
)

# Phase 3: alternative config — unchanged
phase3_fp_alt = Pulse(
    pulse_type="FP_ALT",
    nb_pulses=5,
    ramp_up=400,
    steady_state=600,
    ramp_down=400,
    waiting=3600,
    tritium_fraction=0.5,
    heat_scaling=1.0,
    flux_scaling=1.0,
)

# Baking — unchanged
bake = Pulse(
    pulse_type="BAKE",
    nb_pulses=1,
    ramp_up=151200,
    steady_state=345600,
    ramp_down=108000,
    waiting=10,
    tritium_fraction=0.0,
)

# ── Scenario ──────────────────────────────────────────────────────────────────
scenario = Scenario(
    pulses=[phase1_fp, phase2_fp, phase3_fp_alt, bake],
    baking_temp=493.15,
)
scenario.plasma_data_handling = plasma_data_handling
