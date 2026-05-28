# 03 — Materials

**Previous:** [02 — Input Table](02_input_table.md) | **Next:** [04 — Scenario](04_scenario.md)

---

`FW_example/materials.csv` supplies the transport and trapping parameters for
every material name referenced in `input_table.csv`.

Full parameter definitions: [wiki/Input-Files-Reference.md](../../wiki/Input-Files-Reference.md).

## Tutorial materials

Two materials are defined: **W** (standard) and **W_hightraps** (sensitivity variant with 5× trap density, used in [08_sensitivity_scan.md](08_sensitivity_scan.md)).

### W — standard tungsten

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `Mat_density` | 6.3382 × 10²⁸ m⁻³ | Tungsten atom density |
| `D0` | 2.06 × 10⁻⁷ m²/s | Diffusivity pre-exponential |
| `E_D` | 0.28 eV | Diffusion activation energy |
| `K_R` | 7.94 × 10⁻¹⁷ m⁴/s | Recombination rate coefficient |
| `E_R` | −2 eV | Recombination activation energy |
| `N_traps` | 2 | Number of trap types |

| Trap | `Trap_density` (at. frac.) | `E_k` (eV) | `E_p` (eV) |
|------|--------------------------|-----------|-----------|
| 1 | 1 × 10⁻⁴ | 0.28 | 0.85 |
| 2 | 1 × 10⁻⁴ | 0.28 | 1.00 |

### W_hightraps — sensitivity variant

Same diffusion/recombination parameters as **W**, but trap densities increased
to 5 × 10⁻⁴ (5× higher). This exaggerates retention to probe model sensitivity —
see [08_sensitivity_scan.md](08_sensitivity_scan.md).

## File format reminder

`materials.csv` uses a **horizontal layout**: each material occupies two columns
(key, value), with blank columns between materials. HISP reads as many materials
as it finds column pairs.

```
Material_name, W,,Material_name, W_hightraps
Mat_density, 6.3382E+028,,Mat_density, 6.3382E+028
...
```

Do **not** add a header row or change the key names — they are matched exactly
by the loader.

---

**Next:** [04 — Scenario](04_scenario.md)
