"""
temperature_models.py — Tutorial: First Wall
============================================
Optional custom temperature model for W bins.

HISP ships a built-in W temperature model (Holzner thermal conductivity +
Cu substrate). This file shows how to override it with a simplified
linear profile calibrated to the tutorial geometry.

Only materials listed in TEMPERATURE_MODELS are overridden; all others
continue to use HISP defaults.

See wiki/Temperature-Models.md for full model documentation.
"""

import numpy as np


# ── Calibration data: (heat_load_W/m², T_surface_K, T_rear_K) ────────────────
# Estimated from steady-state 1-D conduction:
#   W: k ≈ 160 W/(m·K), thickness = 6 mm → ΔT = q·L/k
#   Cu substrate: k ≈ 390 W/(m·K), thickness = 2 mm
#   Coolant: 70 °C (343 K)
#
# Example: q = 0.5 MW/m²
#   ΔT_W  = 0.5e6 × 0.006 / 160 = 18.75 K
#   ΔT_Cu = 0.5e6 × 0.002 / 390 =  2.56 K
#   T_rear ≈ 343 + 2.56 ≈ 346 K
#   T_surf ≈ 346 + 18.75 ≈ 365 K
W_HEAT_TEMP_TRIADS = [
    (0.0,    343.0,  343.0),   # no heat load → coolant temperature
    (2.0e5,  369.0,  344.0),   # 0.2 MW/m²
    (5.0e5,  399.0,  346.0),   # 0.5 MW/m²
    (1.0e6,  456.0,  349.0),   # 1.0 MW/m² (FW design limit)
]


def _interpolate_temperatures(heat_flux, triads):
    data = np.asarray(triads, dtype=float)
    order = np.argsort(data[:, 0])
    q     = data[order, 0]
    t_s   = data[order, 1]
    t_r   = data[order, 2]
    return float(np.interp(float(heat_flux), q, t_s)), \
           float(np.interp(float(heat_flux), q, t_r))


def w_temperature_model(x, heat_flux, thickness, **kwargs):
    """
    Linear temperature profile through the W layer.

    Parameters
    ----------
    x         : position(s) in the W layer [m], 0 = plasma-facing surface
    heat_flux  : incident heat load [W/m²]
    thickness  : W layer thickness [m]

    Returns
    -------
    Temperature [K] at position x (scalar or array matching x)
    """
    t_surface, t_rear = _interpolate_temperatures(heat_flux, W_HEAT_TEMP_TRIADS)
    x_arr = np.asarray(x, dtype=float)
    if thickness <= 0:
        return np.full_like(x_arr, t_surface)
    return t_surface + (t_rear - t_surface) * (x_arr / float(thickness))


# ── Export override dict ──────────────────────────────────────────────────────
# Keys must match material names in materials.csv (case-sensitive).
TEMPERATURE_MODELS = {
    "W":          w_temperature_model,
    "W_hightraps": w_temperature_model,  # same geometry, same thermal model
}
