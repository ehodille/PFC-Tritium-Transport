#!/usr/bin/env python3
"""
PFC-Tritium-Transport — interactive local runner.

Prompts for an input folder and a simulation ID, then runs that bin locally
using the same logic as run_on_cluster/run_bin_from_folder.py.

Usage:
    python main.py

You will be asked for:
    - Input folder  : path to the folder containing input_table.csv, scenario.py, etc.
    - Simulation ID : the integer ID of the bin to run (from the 'Sim. ID' column
                      in input_table.csv, or 1-based row number if the column is absent).
"""

import os
import sys

# ---------------------------------------------------------------------------
# Path bootstrap
# ---------------------------------------------------------------------------
_REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

_HISP_SRC = os.path.join(_REPO_ROOT, "hisp", "src")
if _HISP_SRC not in sys.path:
    sys.path.insert(0, _HISP_SRC)

# ---------------------------------------------------------------------------
# Interactive prompts
# ---------------------------------------------------------------------------
print("=" * 60)
print("  PFC-Tritium-Transport — local simulation runner")
print("=" * 60)

input_folder = input("\nInput folder name: ").strip()
if not input_folder:
    print("Error: input folder cannot be empty.")
    sys.exit(1)

sim_id_str = input("Simulation ID    : ").strip()
try:
    sim_id = int(sim_id_str)
except ValueError:
    print(f"Error: simulation ID must be an integer, got '{sim_id_str}'.")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Delegate to run_bin_from_folder.py
# ---------------------------------------------------------------------------
_RUN_ON_CLUSTER = os.path.join(_REPO_ROOT, "run_on_cluster")
if _RUN_ON_CLUSTER not in sys.path:
    sys.path.insert(0, _RUN_ON_CLUSTER)

sys.argv = [
    "run_on_cluster/run_bin_from_folder.py",
    input_folder,
    str(sim_id),
]

runner = os.path.join(_RUN_ON_CLUSTER, "run_bin_from_folder.py")

import runpy
runpy.run_path(runner, run_name="__main__")

