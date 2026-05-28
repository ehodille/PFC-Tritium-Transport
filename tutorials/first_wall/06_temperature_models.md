# 06 — Temperature Models

**Previous:** [05 — Mesh](05_mesh.md) | **Next:** [07 — Running](07_running.md)

---

`FW_example/temperature_models.py` optionally overrides HISP's built-in
temperature models with custom functions calibrated to the tutorial geometry.

Full model documentation: [wiki/Temperature-Models.md](../../wiki/Temperature-Models.md).

## When to use a custom model

HISP has built-in temperature models for W (with/without Cu substrate),
B, and SS. Use a custom model when:
- You have experimental or FEM-derived surface temperature data.
- You want to test sensitivity to a different thermal conductivity assumption.
- You are using a material not in the built-in set.

For this tutorial, the custom W model uses a **simple linear profile** based on
1-D steady-state conduction, calibrated to the 6 mm W + 2 mm Cu + 70 °C
coolant geometry.

## Calibration triads

```python
W_HEAT_TEMP_TRIADS = [
    (0.0,    343.0,  343.0),   # no heat → coolant temperature everywhere
    (2.0e5,  369.0,  344.0),   # 0.2 MW/m²
    (5.0e5,  399.0,  346.0),   # 0.5 MW/m²
    (1.0e6,  456.0,  349.0),   # 1.0 MW/m² (FW design limit)
]
```

Each triad is `(heat_load [W/m²], T_surface [K], T_rear [K])`. HISP
interpolates between triads at each time step.

*Derivation (0.5 MW/m² example):*
$$\Delta T_W = \frac{q \cdot L_W}{k_W} = \frac{5\times10^5 \times 0.006}{160} \approx 19 \text{ K}$$

## The model function interface

```python
def w_temperature_model(x, heat_flux, thickness, **kwargs):
    """
    x         : position in W layer [m], 0 = plasma-facing surface
    heat_flux  : incident heat load [W/m²]
    thickness  : W layer thickness [m]
    **kwargs   : reserved for future HISP arguments
    Returns   : temperature [K] at x (scalar or array)
    """
    t_surface, t_rear = _interpolate_temperatures(heat_flux, W_HEAT_TEMP_TRIADS)
    x_arr = np.asarray(x, dtype=float)
    return t_surface + (t_rear - t_surface) * (x_arr / float(thickness))
```

## Registration

```python
TEMPERATURE_MODELS = {
    "W":           w_temperature_model,
    "W_hightraps": w_temperature_model,  # same geometry, reuse model
}
```

Keys must exactly match material names in `materials.csv`. Any material not
listed here uses HISP's default model.

---

**Next:** [07 — Running](07_running.md)
