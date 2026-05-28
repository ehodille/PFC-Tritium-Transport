# 04 — Scenario

**Previous:** [03 — Materials](03_materials.md) | **Next:** [05 — Mesh](05_mesh.md)

---

`FW_example/scenario_reference.py` defines what happens to the wall over time:
which pulses fire in what order, at what power, with which plasma data.

Detailed `Scenario`/`Pulse` API: [wiki/Pulses-and-Scenarios.md](../../wiki/Pulses-and-Scenarios.md).
Scaling factor mechanics: [wiki/Pulse-Scaling-Factors.md](../../wiki/Pulse-Scaling-Factors.md).

## Structure of the reference scenario

```
phase1_fp  (10 × FP, 50 % power)   →  fp_pulse.dat  × flux_scaling=0.5
phase2_fp  (20 × FP, 100 % power)  →  fp_pulse.dat  × flux_scaling=1.0
phase3_fp_alt (5 × FP_ALT)         →  alternative_fp_pulse.dat
bake       (1 × BAKE, 220 °C)
```

Total time ≈ 35 pulses × (400 + 600 + 400 + 3600) s + baking duration.

## Why two plasma data files?

```python
plasma_data_handling = PlasmaDataHandling(
    pulse_type_to_data={
        "FP":     pd.read_csv(_data_dir / "fp_pulse.dat"),
        "FP_ALT": pd.read_csv(_data_dir / "alternative_fp_pulse.dat"),
    },
)
```

HISP looks up the pulse type string (`"FP"`, `"FP_ALT"`, …) as a key in
`pulse_type_to_data` to retrieve the flux table for that pulse. This means:
- Different pulse types can use **different** spatial flux distributions.
- A single scenario file can represent distinct plasma configurations without
  duplicating the rest of the input.

## Annotated excerpt

```python
phase1_fp = Pulse(
    pulse_type="FP",
    nb_pulses=10,
    ramp_up=400,           # seconds — plasma current ramp
    steady_state=600,      # seconds — flat-top (actual burning phase)
    ramp_down=400,
    waiting=3600,          # ~1 h dwell between pulses
    tritium_fraction=0.5,  # fraction of hydrogen that is tritium
    heat_scaling=0.5,      # multiply heat_total in plasma data by 0.5
    flux_scaling=0.5,      # multiply all particle fluxes by 0.5
)
```

```python
bake = Pulse(
    pulse_type="BAKE",
    nb_pulses=1,
    ramp_up=151200,        # 42 h at 5 °C/h → from 70 °C to 220 °C
    steady_state=345600,   # 4 days at 220 °C
    ramp_down=108000,      # 30 h at 7 °C/h → back to 70 °C
    waiting=10,
    tritium_fraction=0.0,  # no plasma during baking
)

scenario = Scenario(
    pulses=[phase1_fp, phase2_fp, phase3_fp_alt, bake],
    baking_temp=493.15,    # K — required when any BAKE pulse is present
)
scenario.plasma_data_handling = plasma_data_handling
```

> `baking_temp` must be set whenever a `"BAKE"` or `"Bake+GDC"` pulse is in
> the pulse list; HISP raises an error otherwise.

## Paths

`Path(__file__).resolve().parent.parent / "plasma_data"` resolves to the
`plasma_data/` folder regardless of where the runner is called from. This is
the recommended pattern for tutorial/example scenarios.

---

**Next:** [05 — Mesh](05_mesh.md)
