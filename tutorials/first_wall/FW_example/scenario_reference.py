"""
scenario_reference.py — Tutorial: First Wall (reference case)
==============================================================
Standard ITER first-wall scenario: 30 FP pulses across two operation phases
followed by a baking event.

Run with:
    bash run_on_cluster/slurm_folder_jobs.sh tutorials/first_wall/FW_example

See 04_scenario.md for a line-by-line explanation.
"""

from pathlib import Path
from scenario import Scenario, Pulse
import pandas as pd
from plasma_data_handling import PlasmaDataHandling

# ── Plasma data ───────────────────────────────────────────────────────────────
# Paths are relative to this file, so they work regardless of the working
# directory from which the runner is called.
_data_dir = Path(__file__).resolve().parent.parent / "plasma_data"

plasma_data_handling = PlasmaDataHandling(
    pulse_type_to_data={
        # Standard DT H-mode full-power pulses
        "FP":     pd.read_csv(_data_dir / "fp_pulse.dat"),
        # Alternative magnetic configuration (different flux distribution)
        "FP_ALT": pd.read_csv(_data_dir / "alternative_fp_pulse.dat"),
    },
)

# ── Pulses ────────────────────────────────────────────────────────────────────
# Phase 1: ramp-up campaign — reduced power (50 %) for commissioning
phase1_fp = Pulse(
    pulse_type="FP",
    nb_pulses=10,
    ramp_up=400,           # 400 s plasma current ramp
    steady_state=600,      # 600 s flat-top (ITER Q=10 design pulse)
    ramp_down=400,
    waiting=3600,          # ~1 h between pulses
    tritium_fraction=0.5,  # 50/50 DT mix
    heat_scaling=0.5,      # 50 % of reference heat load
    flux_scaling=0.5,      # 50 % of reference particle flux
)

# Phase 2: full-power campaign — standard plasma
phase2_fp = Pulse(
    pulse_type="FP",
    nb_pulses=20,
    ramp_up=400,
    steady_state=600,
    ramp_down=400,
    waiting=3600,
    tritium_fraction=0.5,
    heat_scaling=1.0,
    flux_scaling=1.0,
)

# Phase 3: alternative magnetic config — different flux spatial distribution
#   Uses "FP_ALT" key → reads from alternative_fp_pulse.dat
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

# Baking event — wall heats to baking_temp, no particle flux
bake = Pulse(
    pulse_type="BAKE",
    nb_pulses=1,
    ramp_up=151200,    # 42 h ramp at 5 °C/h
    steady_state=345600,
    ramp_down=108000,  # 30 h cooldown at 7 °C/h
    waiting=10,
    tritium_fraction=0.0,
)

# ── Scenario ──────────────────────────────────────────────────────────────────
scenario = Scenario(
    pulses=[phase1_fp, phase2_fp, phase3_fp_alt, bake],
    baking_temp=493.15,    # 220 °C ITER first-wall baking temperature
)
scenario.plasma_data_handling = plasma_data_handling
