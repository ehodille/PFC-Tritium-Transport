# 02 — Input Table

**Previous:** [01 — Plasma Data](01_plasma_data.md) | **Next:** [03 — Materials](03_materials.md)

---

`FW_example/input_table.csv` defines every simulation bin: geometry, material,
boundary conditions, and solver tolerances.

Full column reference: [wiki/Input-Files-Reference.md](../../wiki/Input-Files-Reference.md).

## Tutorial bins at a glance

| Sim. ID | Bin number | mode | Material | `f` | `ion_scaling_factor`* | BC (plasma-facing) |
|---------|-----------|------|----------|-----|----------------------|-------------------|
| 1 | 0 | wetted | W | 0.88 | 10.4 | Robin (surf. rec. + implant.) |
| 2 | 1 | wetted | W | 0.95 | 7.1 | Dirichlet (analytical approx.) |
| 3 | 2 | shadowed | W | 0.00 | 0.0 | Robin (surf. rec. + implant.) |

\*`ion_scaling_factor = f × S_parent / S_bin`.
For bin 1: 0.88 × 26.0 / 2.2 = 10.4.
The shadowed bin has f = 0 → no particle flux regardless of plasma data.

## Design choices illustrated

**Two different plasma-facing BCs in the same table:**
- Bin 0 uses *Robin* (physically rigorous, solves surface recombination + implantation together).
- Bin 1 uses *Dirichlet — analytical implantation approximation* (faster; applies an analytical expression for the implanted concentration as a Dirichlet condition). `Calculate Implantation Parameters = Yes` in both cases so HISP computes the range and reflection from the energy/angle data.

**Shadowed bin (bin 2):**
- `f = 0.0` zeroes the particle flux; only heat load and diffusion matter.
- `Calculate Implantation Parameters = No` — no flux data to compute from.
- Uses `Dirichlet - 0 concentration` at the rear surface, making it a simple permeation test case alongside the retention study.

## Key column values

```
Thickness (m)    = 0.006   → 6 mm W armour
Cu thickness (m) = 0.002   → 2 mm Cu heat sink (used by thermal model)
rtol             = 1E-12   → tight relative tolerance
atol             = 1E9     → absolute tolerance in atoms/m³
FP max. stepsize = 1000 s  → maximum solver step during plasma pulses
Max. stepsize no FP = 10000 s → maximum step outside plasma pulses
location         = FW      → used by post-processing scripts
```

> **Note on `Sim. ID`:** If omitted, the runner assigns 1-based row numbers.
> Explicit Sim. IDs let you add rows later without renumbering existing checkpoints.

## Rear surface BC and permeation studies

The `BC rear surface` column controls what happens at the back of the material slab.
The two most common choices are:

| BC option | Physical meaning | Use case |
|-----------|-----------------|----------|
| `Neumann - no flux` | Zero permeation flux — tritium cannot leave through the rear | Default for armour tiles backed by Cu/CuCrZr (coolant barrier); gives **inventory** as the output of interest |
| `Dirichlet - 0 concentration` | Concentration fixed at zero at the rear — tritium that reaches the back is immediately removed | Models a perfectly absorbing substrate or a vacuum side; gives **inventory** as well as **permeation flux** as the output of interest |
| `Robin - Surface Recombination` | Recombination-limited release at the rear surface | Models a metallic rear surface where molecules must recombine before release; intermediate between the two above, also valid for both **inventory** and **permeation** studies|

For a permeation study (e.g. measuring tritium breakthrough into a coolant channel):

```
BC Plasma Facing Surface = Robin - Surf. Rec. + Implantation ; Dirichlet - Analyttical implantation approximation
BC rear surface          = Dirichlet - 0 concentration ; Robin - Surface Recombination
```

The time-dependent permeation flux at the rear is then available in the output JSON
as the gradient of the concentration at x = L. Compare this to the implanted flux
to quantify the permeation fraction.

> In the tutorial `FW_example`, bins 0 and 1 (wetted) use `Neumann - no flux` at the rear,
> so retention is the output of interest. Bin 2 (shadowed, `f=0`) uses `Dirichlet - 0 concentration`
> at the rear — a low-cost permeation test case since the negligible incoming flux keeps
> runtimes short while demonstrating the BC.

---

**Next:** [03 — Materials](03_materials.md)
