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

```python
graded_vertices(L=bin.thickness, h0=5e-10, r=1.03)
```

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `L` | 6 mm (from `input_table.csv`) | full W layer |
| `h0` | 5 × 10⁻¹⁰ m (0.5 nm) | resolves the ~3 nm implantation range |
| `r` | 1.03 | gentle growth; ~130 nodes across 6 mm |

If a bin's `sim_id` is not in `BINS_MESHES`, HISP generates a default mesh
with `h0=1e-10, r=1.05`. The explicit definition here gives slightly coarser
but sufficient resolution with faster runtime.

## Key snippet from `mesh.py`

```python
for _bin in _reactor.bins:
    _mesh = graded_vertices(L=_bin.thickness, h0=5e-10, r=1.03)
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
