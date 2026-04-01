"""Optional custom temperature models for PFC-Tritium-Transport.

This file is discovered automatically by the folder runner.
Only listed materials are overridden; all others use HISP defaults.
"""

import numpy as np

# Triads: (heat_load_Wm2, surface_temperature_K, rear_temperature_K)
# Replace these example values with your SS calibration data.
SS_HEAT_TEMP_TRIADS = [
    (0.0, 343.0, 343.0),
    (3.5e5, 523.15, 448.15),
    (1.0e6, 700.0, 600.0),
]


def _interp_surface_and_rear(heat_flux, triads):
    data = np.asarray(triads, dtype=float)
    order = np.argsort(data[:, 0])
    q = data[order, 0]
    t_surface = data[order, 1]
    t_rear = data[order, 2]

    ts = float(np.interp(float(heat_flux), q, t_surface))
    tr = float(np.interp(float(heat_flux), q, t_rear))
    return ts, tr


def ss_temperature_model(x, heat_flux, thickness, **kwargs):
    """Custom SS temperature model.

    Inputs are provided by HISP at each evaluation time:
    - `heat_flux` comes from scenario + plasma_data_handling
    - `x` is local position in material (scalar or array)
    - `thickness` is bin thickness

    Per-bin behavior:
    - If `bin.sim_id` matches a key in a per-sim dict, use that sim-specific
      calibration.
    - Otherwise, use `SS_HEAT_TEMP_TRIADS` as default.

    Model:
    1) Interpolate surface and rear temperatures from triads vs heat load.
    2) Build linear profile inside material between surface and rear.
    """
    t_surface, t_rear = _interp_surface_and_rear(heat_flux, SS_HEAT_TEMP_TRIADS)

    x_arr = np.asarray(x, dtype=float)
    if thickness <= 0:
        return np.full_like(x_arr, t_surface)

    return t_surface + (t_rear - t_surface) * (x_arr / float(thickness))


TEMPERATURE_MODELS = {
    "SS": ss_temperature_model,
}
