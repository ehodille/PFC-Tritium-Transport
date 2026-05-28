# 05 — Mesh

**Previous:** [04 — Scenario](04_scenario.md) | **Next:** [06 — Temperature Models](06_temperature_models.md)

---

`FW_example/mesh.py` builds a graded 1-D mesh for each bin. Mesh resolution
directly affects accuracy near the plasma-facing surface (implantation zone)
and simulation runtime.

Full mesh documentation: [wiki/Input-Files-Reference.md](../../wiki/Input-Files-Reference.md#meshpy).

## Graded mesh concept

```
Surface (x=0)                                     Rear (x=L)
  |—·—·——·————·———————·————————————————————————————|
  fine ←————————————————————————— coarse ——————————→
```

The first element has size `h0`; each subsequent element is `r` times larger.

## Tutorial mesh parameters

The tutorial uses two mesh strategies depending on the plasma-facing BC:

| BC (plasma face) | Mesh type | `h0` | `r` | Rationale |
|-----------------|-----------|------|-----|-----------|
| Dirichlet (analytical) | one-sided graded | 1 × 10⁻⁸ m | 1.03 | concentration prescribed at x=0; only surface resolution matters |
| Robin (surf. rec.) | symmetric graded | 1 × 10⁻⁸ m | 1.03 | flux-based BC active at both surfaces; refine at x=0 and x=L |

## Key snippet from `mesh.py`

```python
for _bin in _reactor.bins:
    _bc_plasma = _bin.bin_configuration.bc_plasma_facing_surface
    if "Dirichlet" in _bc_plasma:
        # Concentration prescribed at x=0 — refine only at the surface
        _mesh = graded_vertices(L=_bin.thickness, h0=1e-8, r=1.03)
    else:
        # Robin BC — flux condition at both surfaces; refine at both ends
        _mesh = symmetric_graded_vertices(L=_bin.thickness, h0=1e-8, r=1.03)
    BINS_MESHES[_bin.sim_id] = MeshBin(sim_id=_bin.sim_id, mesh=_mesh)
```

> `_reactor` is loaded from `input_table.csv` at import time using the
> `INPUT_DIR_CONTEXT` environment variable (set by the SLURM runner) or, when
> absent, the directory containing `mesh.py` itself.

## Choosing h0 and r

- **h0**: should be ≤ implantation range / 3. For W at typical FW energies
  (85–165 eV), range ≈ 2–5 nm → h0 ≤ 1 nm is safe. Smaller h0 increases
  node count and runtime.
- **r**: values 1.02–1.05 are typical. Lower → more nodes, higher accuracy near
  rear; higher → faster but coarser rear resolution.
- For shadowed bins (no implantation), a coarser mesh (larger h0, larger r) is
  acceptable.

---

**Next:** [06 — Temperature Models](06_temperature_models.md)
