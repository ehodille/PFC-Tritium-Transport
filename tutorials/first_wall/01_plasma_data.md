# 01 — Plasma Data Files

**Previous:** [README](README.md) | **Next:** [02 — Input Table](02_input_table.md)

---

The plasma data files provide per-bin particle fluxes, implantation parameters,
and heat loads for each pulse type. Each row corresponds to one bin (identified
by `Bin_Index`, which equals the `Bin number` in `input_table.csv`).

## Format

```
Bin_Index, Flux_Ion, Flux_Atom, E_ion, E_atom, alpha_ion, alpha_atom, heat_total, heat_ion
```

| Column | Units | Description |
|--------|-------|-------------|
| `Bin_Index` | — | Integer bin index (must match `Bin number` in input_table) |
| `Flux_Ion` | m⁻²s⁻¹ | Ion particle flux |
| `Flux_Atom` | m⁻²s⁻¹ | Neutral atom flux (charge-exchange) |
| `E_ion` | eV | Mean ion impact energy |
| `E_atom` | eV | Mean atom impact energy |
| `alpha_ion` | degrees | Mean ion incidence angle (from surface normal) |
| `alpha_atom` | degrees | Mean atom incidence angle |
| `heat_total` | W/m² | Total heat load |
| `heat_ion` | W/m² | Ion-carried heat load component |

> The `Bin_Index` value is also used by `ImplantationCalculator` to look up
> implantation range and reflection coefficients. Both column name and row
> position must match.

## Tutorial data

This tutorial uses **three FW bins** (indices 0–2):

| Bin | Location | `Flux_Ion` (m⁻²s⁻¹) | `heat_total` (W/m²) | Notes |
|-----|----------|----------------------|----------------------|-------|
| 0 | mid-FW | 2.5 × 10¹⁹ | 8.5 × 10⁴ | moderately loaded |
| 1 | upper FW | 6.8 × 10¹⁹ | 1.4 × 10⁵ | highest flux |
| 2 | lower FW | 5.0 × 10¹⁶ | 1.2 × 10⁴ | shadowed, negligible flux |

Two files are provided:

- [`plasma_data/fp_pulse.dat`](plasma_data/fp_pulse.dat) — standard DT H-mode
- [`plasma_data/alternative_fp_pulse.dat`](plasma_data/alternative_fp_pulse.dat) — alternative magnetic configuration (shifted peak flux, ~30 % higher on bin 0)

The scenario maps each file to a pulse type string (`"FP"` and `"FP_ALT"`
respectively) so HISP selects the right flux table during each pulse phase.

## How these values were chosen

ITER first-wall design parameters (for reference):
- Peak particle flux: ~10²⁰ m⁻²s⁻¹ at the equatorial FW
- Design heat load limit: 1 MW/m²
- Typical ion energy at FW: 50–200 eV (sheath potential ~3T_e)
- Typical incidence angle: 60–80° (oblique relative to surface)

These tutorial values are representative but simplified; the actual SOLPS-ITER
binned data for real simulations is in `data/` and `binning/`.

---

**Next:** [02 — Input Table](02_input_table.md)
