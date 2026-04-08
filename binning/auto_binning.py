#!/usr/bin/env python3
"""
Automatic binning programme.

Reads SOLPS output files (segment coordinates + plasma data) and reactor
wall coordinate files, creates bins, maps SOLPS values onto them, and
produces output data files and diagnostic plots.

Two kinds of reactor coordinate files are supported:
  - **fixed_bins/**  : Pre-defined bins (e.g. First Wall).  SOLPS values
        are mapped directly onto these bins (length-weighted average).
        No splitting, shadow-filtering or merging is applied.
  - **coordinates_to_be_binned/**  : Raw wall segments (e.g. divertor).
        These are processed: shadow-filtered, subdivided into fine bins,
        SOLPS-mapped, and then merged back with a relative tolerance.

Directory layout expected
─────────────────────────
    binning/
    ├── auto_binning.py          (this script)
    ├── SOLPS_data/
    │   ├── wlld_full-shd.dat   (SOLPS plasma data, 37 columns)
    │   └── wlly_full-shd.dat   (SOLPS segment coordinates)
    ├── Reactor_coordinates/
    │   ├── fixed_bins/          (one .dat file per region, FW format)
    │   └── coordinates_to_be_binned/   (one .dat file per region, polyline)
    └── output/                  (all generated files & plots)

Usage:
    python auto_binning.py
"""

import os
import sys
import csv
import glob as _glob
import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict

# ──────────────────────── Configuration ──────────────────────────────────────

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# ── Input directories ──
SOLPS_DIR       = os.path.join(BASE_DIR, "SOLPS_data")
REACTOR_DIR     = os.path.join(BASE_DIR, "Reactor_coordinates")
FIXED_BINS_DIR  = os.path.join(REACTOR_DIR, "fixed_bins")
COORDS_TO_BIN_DIR = os.path.join(REACTOR_DIR, "coordinates_to_be_binned")

# ── Output directory ──
OUTPUT_DIR = os.path.join(BASE_DIR, "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── SOLPS input files ──
# Auto-detect wlld file (may be wlld_full-shd.dat or wlld_full-shd.NNN.dat)
_wlld_candidates = sorted(_glob.glob(os.path.join(SOLPS_DIR, "wlld_full-shd*.dat")),
                           key=len)
WLLD_PATH = _wlld_candidates[0] if _wlld_candidates else os.path.join(SOLPS_DIR, "wlld_full-shd.dat")

_wlly_candidates = sorted(_glob.glob(os.path.join(SOLPS_DIR, "wlly_full-shd*.dat")),
                           key=len)
WLLY_PATH = _wlly_candidates[0] if _wlly_candidates else os.path.join(SOLPS_DIR, "wlly_full-shd.dat")

# ── Reactor coordinate files (auto-discovered) ──
# fixed_bins: every .dat file in fixed_bins/ is treated as a pre-defined bin set
FIXED_BIN_FILES = sorted(_glob.glob(os.path.join(FIXED_BINS_DIR, "*.dat")))

# coordinates_to_be_binned: every .dat file is a polyline to be processed
COORDS_TO_BIN_FILES = sorted(_glob.glob(os.path.join(COORDS_TO_BIN_DIR, "*.dat")))

# ── Output paths ──
OUT_COORDS         = os.path.join(OUTPUT_DIR, "new_bins_coords.dat")
OUT_DATA           = os.path.join(OUTPUT_DIR, "new_bins_data.dat")
OUT_PLOT           = os.path.join(OUTPUT_DIR, "auto_binning_plot.png")
OUT_PREMRGE_COORDS = os.path.join(OUTPUT_DIR, "pre_merge_bins_coords.dat")
OUT_PREMRGE_DATA   = os.path.join(OUTPUT_DIR, "pre_merge_bins_data.dat")
OUT_FW_COORDS      = os.path.join(OUTPUT_DIR, "mapped_fixed_bins_coords.dat")
OUT_FW_DATA        = os.path.join(OUTPUT_DIR, "fixed_bins_data.dat")
OUT_BINNED_FLUX_CSV = os.path.join(OUTPUT_DIR, "Binned_Flux_Data.dat")

# -- Filtering & binning parameters --

# Maximum distance (m) from a divertor point to the nearest SOLPS segment
# centre for that divertor region to be considered "illuminated" (not shadowed).
MAX_SHADOW_DIST = 0.1  # metres

# Maximum distance (m) for an SOLPS segment to contribute to a bin's
# weighted flux mapping.
MAX_MAPPING_DIST = 0.05  # metres

# Maximum angle (degrees) between an SOLPS segment tangent and the local
# divertor tangent for the SOLPS segment to be considered relevant (filters
# out perpendicular or misaligned SOLPS segments).
MAX_ANGLE_DEG = 10.0

# Minimum bin arc-length (m). Set to 0 to disable min-length enforcement.
MIN_BIN_LENGTH = 0.0  # metres  (disabled)

# Relative tolerance for merging adjacent bins with similar fluxes.
MERGE_RTOL = 0.20  # 20 %

# Tighter tolerance for bins containing sub-bins with heat flux > threshold
MERGE_RTOL_HOT = 0.10   # 10 % (Wtot only)
MERGE_HOT_THRESH = 2e6  # 2 MW/m²

# Maximum allowed bin arc-length after merging (m). Prevents over-merging
# of long similar-flux regions into a single oversized bin.
MAX_BIN_LENGTH = 0.3  # metres

# Minimum edge bin arc-length (m). After merging, if the first or last bin
# is shorter than this, it is absorbed into its neighbor.
MIN_EDGE_BIN_LENGTH = 0.04  # metres (4 cm)

# Number of sub-sample points per divertor segment edge for high-resolution
# mapping from SOLPS to new bins.
SUBSAMPLE_N = 100000

# Column names in wlld_full-shd.dat  (37 columns)
WLLD_COLUMNS = [
    "x", "Wtot", "Wrad", "Wpart", "Wpls", "Wneu",
    "Wheat", "Wpot", "Whtpl", "Wptpl",
    "ne", "Te", "Ti",
    "flxi_D", "Eavi_D", "flxi_He", "Eavi_He", "flxi_Ne", "Eavi_Ne",
    "flxa_D", "Eava_D", "na_D", "pa_D",
    "flxa_He", "Eava_He", "na_He", "pa_He",
    "flxa_Ne", "Eava_Ne", "na_Ne", "pa_Ne",
    "flxm_D2", "Eavm_D2", "nm_D2", "pm_D2",
    "ersne_W", "ersng_W",
]

MERGE_HEAT_KEY = "Wtot"
MERGE_DFLUX_KEYS = ("flxi_D", "flxa_D")

# Absolute tolerance for merging – no longer used but kept for compat
MERGE_ATOL = 0.0  # disabled

# ──────────────────────── Data Structures ────────────────────────────────────

@dataclass
class SOLPSSegment:
    """One SOLPS wall segment with geometry and plasma data."""
    r_start: float
    z_start: float
    r_centre: float
    z_centre: float
    r_end: float
    z_end: float
    x: float          # arc-length coordinate along SOLPS plot zone
    seg_id: int
    data: Dict[str, float] = field(default_factory=dict)


@dataclass
class NewBin:
    """A new divertor bin."""
    bin_id: int
    r_start: float
    z_start: float
    r_end: float
    z_end: float
    arc_length: float  # metres
    data: Dict[str, float] = field(default_factory=dict)
    n_solps_contributors: int = 0


# ──────────────────────── Parsing Functions ──────────────────────────────────

def parse_wlly(path: str) -> List[SOLPSSegment]:
    """Parse wlly_full-shd.dat → list of SOLPSSegment (geometry only)."""
    segments = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                r = float(parts[0])
                z = float(parts[1])
                rc = float(parts[2])
                zc = float(parts[3])
                x = float(parts[4])
                seg_id = int(parts[5]) if len(parts) >= 6 else -1
            except (ValueError, IndexError):
                continue
            # Endpoint by reflection about centre: end = 2*centre - start
            r_end = 2.0 * rc - r
            z_end = 2.0 * zc - z
            segments.append(SOLPSSegment(
                r_start=r, z_start=z,
                r_centre=rc, z_centre=zc,
                r_end=r_end, z_end=z_end,
                x=x, seg_id=seg_id,
            ))
    return segments


def parse_wlld(path: str) -> Dict[float, Dict[str, float]]:
    """Parse wlld_full-shd.dat → dict keyed by x coordinate."""
    data = {}
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 37:
                continue
            try:
                vals = [float(p) for p in parts[:37]]
            except ValueError:
                continue
            row = {col: vals[i] for i, col in enumerate(WLLD_COLUMNS)}
            data[row["x"]] = row
    return data

def parse_wlld_rows(path: str) -> List[Dict[str, float]]:
    """
    Parse wlld_full-shd.dat and return a list of dicts preserving file row order.
    Each dict contains the 37 WLLD_COLUMNS.
    """
    rows: List[Dict[str, float]] = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 37:
                continue
            try:
                vals = [float(p) for p in parts[:37]]
            except ValueError:
                continue
            row = {col: vals[i] for i, col in enumerate(WLLD_COLUMNS)}
            rows.append(row)
    return rows


def parse_divertor(path: str) -> List[List[Tuple[float, float]]]:
    """Parse divertor.dat → list of polylines (each a list of (R, Z) points).
    Groups are separated by blank lines."""
    groups: List[List[Tuple[float, float]]] = []
    current: List[Tuple[float, float]] = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                if current:
                    groups.append(current)
                    current = []
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                r = float(parts[0])
                z = float(parts[1])
                current.append((r, z))
            except ValueError:
                continue
    if current:
        groups.append(current)
    return groups


# ──────────────────────── Geometry Helpers ────────────────────────────────────

def seg_length(r1: float, z1: float, r2: float, z2: float) -> float:
    return np.hypot(r2 - r1, z2 - z1)


def seg_tangent(r1: float, z1: float, r2: float, z2: float) -> np.ndarray:
    """Unit tangent vector of segment (r1,z1)->(r2,z2)."""
    d = np.array([r2 - r1, z2 - z1], dtype=float)
    n = np.linalg.norm(d)
    if n < 1e-15:
        return np.array([1.0, 0.0])
    return d / n
def angle_between_segments(t1: np.ndarray, t2: np.ndarray) -> float:
    """Angle (degrees) between two tangent vectors (unsigned, 0-90)."""
    cos_a = np.clip(abs(np.dot(t1, t2)), 0.0, 1.0)
    return np.degrees(np.arccos(cos_a))
def has_clear_los_to_any_solps(
    mid_r: float, mid_z: float,
    solps_centres: np.ndarray,        # shape (N, 2)
    div_edges: List[Tuple[float, float, float, float]],
    search_radius: float,
    eps: float = 1e-6
) -> bool:
    """
    True if there exists at least one SOLPS center within search_radius such that
    the straight segment Centre->Mid does not properly intersect any divertor edge
    (i.e., clear line-of-sight).

    We ignore intersections within 'eps' of either endpoint to prevent false
    positives from touching the wall exactly at P or Q.
    """
    if solps_centres.size == 0:
        return False

    P = np.array([mid_r, mid_z], dtype=float)
    # shortlist centres by distance
    d2 = np.sum((solps_centres - P[None, :])**2, axis=1)
    mask = d2 <= (search_radius * search_radius)
    if not np.any(mask):
        return False

    candidates = solps_centres[mask]
    for C in candidates:
        cx, cy = float(C[0]), float(C[1])
        # Quick AABB for PQ to prune edges
        minx = min(cx, mid_r) - eps; maxx = max(cx, mid_r) + eps
        miny = min(cy, mid_z) - eps; maxy = max(cy, mid_z) + eps

        blocked = False
        for (ax, ay, bx, by) in div_edges:
            # AABB overlap test
            if (max(min(ax, bx), minx) > min(max(ax, bx), maxx)) or \
               (max(min(ay, by), miny) > min(max(ay, by), maxy)):
                continue
            if _segments_intersect_proper(cx, cy, mid_r, mid_z, ax, ay, bx, by, eps=eps):
                blocked = True
                break

        if not blocked:
            return True  # found a visible centre

    return False


# ──────────────────────── New: Self-shadow occlusion at sub-bin level ───────

def collect_divertor_edges(div_groups: List[List[Tuple[float, float]]]
                           ) -> List[Tuple[float, float, float, float]]:
    """Flatten divertor polylines into a list of edges (ax, ay, bx, by)."""
    edges = []
    for grp in div_groups:
        if len(grp) < 2:
            continue
        for i in range(len(grp) - 1):
            (ax, ay) = grp[i]
            (bx, by) = grp[i + 1]
            edges.append((ax, ay, bx, by))
    return edges


def _segments_intersect_proper(px: float, py: float, qx: float, qy: float,
                               ax: float, ay: float, bx: float, by: float,
                               eps: float = 1e-9) -> bool:
    """
    Proper intersection between P->Q and A->B:
    True if they intersect strictly inside both segments (ignores contacts
    within 'eps' of endpoints).
    """
    rx, ry = (qx - px), (qy - py)         # PQ
    sx, sy = (bx - ax), (by - ay)         # AB
    denom = rx * sy - ry * sx
    apx, apy = (ax - px), (ay - py)

    if abs(denom) < eps:
        # Parallel/colinear → treat as not blocking (no proper crossing).
        # (If you want colinear-overlap to block, add extra checks here.)
        return False

    # Solve P + t r = A + u s
    t = (apx * sy - apy * sx) / denom
    u = (apx * ry - apy * rx) / denom

    return (eps < t < 1.0 - eps) and (eps < u < 1.0 - eps)


def _closest_point_on_segment(px: float, py: float,
                              ax: float, ay: float,
                              bx: float, by: float) -> Tuple[float, float, float]:
    """
    Closest point Q on segment A->B to point P. Returns (qx, qy, t), with t in [0,1].
    """
    vx, vy = (bx - ax), (by - ay)
    wx, wy = (px - ax), (py - ay)
    vv = vx*vx + vy*vy
    if vv <= 1e-30:
        return ax, ay, 0.0
    t = (wx*vx + wy*vy) / vv
    t = max(0.0, min(1.0, t))
    qx, qy = ax + t * vx, ay + t * vy
    return qx, qy, t


def point_to_segment_dist(px: float, py: float,
                          ax: float, ay: float,
                          bx: float, by: float) -> float:
    """Distance from point P to segment A->B."""
    qx, qy, _ = _closest_point_on_segment(px, py, ax, ay, bx, by)
    return float(np.hypot(px - qx, py - qy))


def densify_polyline(pts: List[Tuple[float, float]], n_per_edge: int
                     ) -> List[Tuple[float, float]]:
    """Return a densified polyline inserting `n_per_edge` equally spaced
    points along each edge (excluding the trailing endpoint except for the
    final segment where the original last point is kept).
    """
    if not pts:
        return []
    if n_per_edge <= 0:
        return pts[:]

    out: List[Tuple[float, float]] = []
    for i in range(len(pts) - 1):
        ax, ay = pts[i]
        bx, by = pts[i + 1]
        for k in range(n_per_edge):
            t = k / float(n_per_edge)
            out.append((ax + t * (bx - ax), ay + t * (by - ay)))
    out.append(pts[-1])
    return out

    from typing import List, Tuple
import numpy as np

def snap_discontinuity_endpoints_to_near_segments(
    bins: List[NewBin],
    *,
    near_tol: float = 0.01,           # 1 cm: max distance P->segment to allow snapping
    connect_eps: float = 1e-6,        # when two endpoints are already “the same”
    endpoint_margin: float = 0.01,    # 1 cm: Q must be at least this far from the other segment’s endpoints (i.e., interior)
    drop_if_shorter_than: float = 1e-7,
    max_iters: int = 6
) -> Tuple[List[NewBin], int, List[Tuple[float, float]]]:
    """
    Find dangling endpoints (not shared with another segment). For each dangling endpoint P,
    find the closest point Q on *another* segment (distance to the segment, not to its endpoints).
    If dist(P,Q) <= near_tol and Q is *interior* to that segment (>= endpoint_margin from its
    endpoints), then:
      - snap P -> Q (shortening its segment),
      - trim the other segment so Q becomes its endpoint (delete the shorter tail of that segment).

    Repeat for a few passes to resolve cascades. Remove degenerate segments.

    Returns
    -------
    bins_out : trimmed bins
    n_snaps  : number of snaps applied
    snap_vertices : list of (R, Z) tuples – the new vertices created by each snap
    """

    def _len(bn: NewBin) -> float:
        return seg_length(bn.r_start, bn.z_start, bn.r_end, bn.z_end)

    def _closest_point_on_segment(px: float, py: float,
                                  ax: float, ay: float,
                                  bx: float, by: float) -> Tuple[float, float, float, float]:
        """
        Closest point Q on segment A->B to point P. Returns (qx, qy, t, d2),
        with t in [0,1] and d2 the squared distance |PQ|^2.
        """
        vx, vy = (bx - ax), (by - ay)
        wx, wy = (px - ax), (py - ay)
        vv = vx*vx + vy*vy
        if vv <= 1e-30:
            # degenerate segment: closest is A
            qx, qy = ax, ay
            return qx, qy, 0.0, (px - qx)**2 + (py - qy)**2
        t = (wx*vx + wy*vy) / vv
        t = max(0.0, min(1.0, t))
        qx, qy = ax + t*vx, ay + t*vy
        d2 = (px - qx)**2 + (py - qy)**2
        return qx, qy, t, d2

    def _endpoint_connected_to_any(x: float, y: float, self_idx: int) -> bool:
        """True if (x,y) coincides (within connect_eps) with ANY other bin endpoint."""
        for j, bj in enumerate(out):
            if j == self_idx:
                continue
            if (seg_length(x, y, bj.r_start, bj.z_start) <= connect_eps or
                seg_length(x, y, bj.r_end,   bj.z_end)   <= connect_eps):
                return True
        return False

    # Work on a copy (we’ll mutate geometries)
    out = [NewBin(
        bin_id=bn.bin_id,
        r_start=bn.r_start, z_start=bn.z_start,
        r_end=bn.r_end,     z_end=bn.z_end,
        arc_length=bn.arc_length,
        data=dict(bn.data),
        n_solps_contributors=bn.n_solps_contributors,
    ) for bn in bins]

    n_snaps_total = 0
    snap_vertices: List[Tuple[float, float]] = []
    tol2 = near_tol * near_tol

    for _ in range(max_iters):
        n = len(out)
        if n == 0:
            break

        changed = False

        # Precompute aabbs for pruning point→segment tests
        aabbs = []
        for b in out:
            x0, y0, x1, y1 = b.r_start, b.z_start, b.r_end, b.z_end
            minx, maxx = (x0, x1) if x0 <= x1 else (x1, x0)
            miny, maxy = (y0, y1) if y0 <= y1 else (y1, y0)
            aabbs.append((minx - near_tol, miny - near_tol, maxx + near_tol, maxy + near_tol))

        # Scan endpoints that are discontinuities
        # We restart the scan after each snap to keep logic simple/stable.
        snapped_this_pass = False

        for i, bi in enumerate(out):
            for is_start in (True, False):
                px = bi.r_start if is_start else bi.r_end
                py = bi.z_start if is_start else bi.z_end

                # Skip if already connected to any other endpoint
                if _endpoint_connected_to_any(px, py, self_idx=i):
                    continue

                # Find closest segment (excluding self) by point→segment distance
                best = None   # (j, qx, qy, d2)
                for j, bj in enumerate(out):
                    if j == i:
                        continue

                    # AABB prune for point-to-segment: point must be near bbox
                    minx, miny, maxx, maxy = aabbs[j]
                    if not (minx <= px <= maxx and miny <= py <= maxy):
                        continue

                    qx, qy, t, d2 = _closest_point_on_segment(px, py,
                                                              bj.r_start, bj.z_start,
                                                              bj.r_end,   bj.z_end)
                    if d2 <= tol2:
                        # Q should be interior to bj (far from its endpoints)
                        dQ_start = seg_length(qx, qy, bj.r_start, bj.z_start)
                        dQ_end   = seg_length(qx, qy, bj.r_end,   bj.z_end)
                        if dQ_start >= endpoint_margin and dQ_end >= endpoint_margin:
                            if best is None or d2 < best[3]:
                                best = (j, qx, qy, d2)

                if best is None:
                    continue

                # Apply snap at best (closest) segment
                j, qx, qy, _ = best

                # 1) Snap our dangling endpoint P -> Q
                if is_start:
                    out[i].r_start = qx; out[i].z_start = qy
                else:
                    out[i].r_end   = qx; out[i].z_end   = qy
                out[i].arc_length = seg_length(out[i].r_start, out[i].z_start,
                                               out[i].r_end,   out[i].z_end)

                # 2) Trim the other segment so Q becomes its endpoint:
                #    delete the *shorter* tail to Q (start or end).
                bj = out[j]
                dQ_start = seg_length(qx, qy, bj.r_start, bj.z_start)
                dQ_end   = seg_length(qx, qy, bj.r_end,   bj.z_end)
                if dQ_start <= dQ_end:
                    out[j].r_start = qx; out[j].z_start = qy
                else:
                    out[j].r_end   = qx; out[j].z_end   = qy
                out[j].arc_length = seg_length(out[j].r_start, out[j].z_start,
                                               out[j].r_end,   out[j].z_end)

                snap_vertices.append((qx, qy))
                n_snaps_total += 1
                snapped_this_pass = True
                changed = True
                break  # break endpoint loop to restart scanning after mutation

            if snapped_this_pass:
                break  # restart outer scan

        # Remove degenerate segments and continue if anything changed
        if changed:
            out = [b for b in out if _len(b) > drop_if_shorter_than]
        else:
            break

    # Final refresh & renumber
    for b in out:
        b.arc_length = seg_length(b.r_start, b.z_start, b.r_end, b.z_end)
    for k, b in enumerate(out, start=1):
        b.bin_id = k

    return out, n_snaps_total, snap_vertices

# ──────────────────────── Step 1: Link SOLPS geometry + data ─────────────────

def link_solps_data(segments: List[SOLPSSegment],
                    data_by_x: Dict[float, Dict[str, float]]) -> None:
    """Attach plasma data from wlld to each SOLPS segment (matched on x).

    Always assigns data from the nearest x-key.  Warns when the gap
    exceeds 0.01 (rare: happens for the few wlly segments that have
    no corresponding wlld row).
    """
    x_keys = np.array(sorted(data_by_x.keys()))
    n_approx = 0
    for seg in segments:
        # Find closest x
        idx = np.argmin(np.abs(x_keys - seg.x))
        gap = abs(x_keys[idx] - seg.x)
        seg.data = data_by_x[x_keys[idx]]
        if gap > 0.01:
            n_approx += 1
    if n_approx:
        print(f"  NOTE: {n_approx} segment(s) matched to nearest wlld x "
              f"with gap > 0.01 (approximate data).")

def link_solps_data_by_row(segments: List[SOLPSSegment],
                           rows: List[Dict[str, float]]) -> None:
    """
    Attach WLLD data to each SOLPS segment by row index:
      segment[i].data = rows[i]
    If counts mismatch, assign for the overlapping range and warn.
    """
    n_seg = len(segments)
    n_row = len(rows)
    n = min(n_seg, n_row)

    for i in range(n):
        segments[i].data = rows[i]

    if n_seg != n_row:
        print(f"WARNING: Row count mismatch: segments={n_seg}, wlld rows={n_row}. "
              f"Linked first {n} rows by index. Unmatched items will have default data.")

        # Optional: initialize unmatched segments to NaN (or 0.0 if you prefer)
        for i in range(n, n_seg):
            segments[i].data = {col: float("nan") for col in WLLD_COLUMNS}



# ──────────────────────── Step 2: Filter SOLPS segments ──────────────────────

def filter_solps_for_divertor(segments: List[SOLPSSegment],
                              div_groups: List[List[Tuple[float, float]]],
                              max_dist: float = MAX_MAPPING_DIST,
                              max_angle: float = MAX_ANGLE_DEG,
                              ) -> List[SOLPSSegment]:
    """
    Keep only SOLPS segments that are:
      1. Within max_dist of the divertor,
      2. Have a tangent roughly aligned with the local divertor tangent (< max_angle).
    """
    # Build dense divertor points + tangents
    div_pts = []
    div_tangents = []
    for grp in div_groups:
        if len(grp) < 2:
            continue
        for i in range(len(grp) - 1):
            r1, z1 = grp[i]
            r2, z2 = grp[i + 1]
            mid = ((r1 + r2) / 2.0, (z1 + z2) / 2.0)
            div_pts.append(mid)
            div_tangents.append(seg_tangent(r1, z1, r2, z2))
    div_pts = np.array(div_pts)
    div_tangents = np.array(div_tangents)

    if len(div_pts) == 0:
        return []

    kept = []
    for seg in segments:
        # Distance from SOLPS segment centre to nearest divertor edge midpoint
        sc = np.array([seg.r_centre, seg.z_centre])
        dists = np.linalg.norm(div_pts - sc, axis=1)
        min_idx = np.argmin(dists)
        min_dist = dists[min_idx]

        if min_dist > max_dist:
            continue

        # Angle check
        solps_tan = seg_tangent(seg.r_start, seg.z_start, seg.r_end, seg.z_end)
        angle = angle_between_segments(solps_tan, div_tangents[min_idx])
        if angle > max_angle:
            continue

        kept.append(seg)

    return kept


# ──────────────────────── Step 3: Create initial bins on divertor ─────────────

def create_initial_bins(div_groups: List[List[Tuple[float, float]]],
                        solps_segments: List[SOLPSSegment],
                        max_shadow_dist: float = MAX_SHADOW_DIST,
                        ) -> List[NewBin]:
    """
    Walk along each divertor polyline edge and keep only those edges whose
    midpoint has a clear line-of-sight to at least one (filtered) SOLPS
    segment centre within max_shadow_dist. Each kept edge becomes one bin.
    """
    if not solps_segments:
        print("WARNING: No filtered SOLPS segments → no bins created.")
        return []

    # Collect SOLPS centres (N,2)
    solps_centres = np.array([[s.r_centre, s.z_centre] for s in solps_segments], dtype=float)

    # Collect all divertor edges for occlusion tests
    div_edges = collect_divertor_edges(div_groups)

    bins: List[NewBin] = []
    bid = 1
    for grp in div_groups:
        if len(grp) < 2:
            continue
        for i in range(len(grp) - 1):
            r1, z1 = grp[i]
            r2, z2 = grp[i + 1]
            mid_r = (r1 + r2) / 2.0
            mid_z = (z1 + z2) / 2.0

            # NEW: visibility/occlusion test to enforce shadowing
            visible = has_clear_los_to_any_solps(
                mid_r, mid_z,
                solps_centres=solps_centres,
                div_edges=div_edges,
                search_radius=max_shadow_dist,
                eps=1e-6,
            )
            if not visible:
                # Shadowed (no clear LoS to any relevant SOLPS centre) → skip
                continue

            length = seg_length(r1, z1, r2, z2)
            bins.append(NewBin(
                bin_id=bid, r_start=r1, z_start=z1,
                r_end=r2, z_end=z2, arc_length=length,
            ))
            bid += 1

    return bins


# ──────────────────────── Step 4: Subdivide long bins ────────────────────────

from typing import List, Optional

def subdivide_bins(
    bins: List[NewBin],
    solps_segments: List[SOLPSSegment],  # unused; kept for signature compatibility
    target_length: Optional[float] = None,
    fixed_vertices: Optional[List[Tuple[float, float]]] = None,
) -> List[NewBin]:
    """
    Subdivide the *entire* sequence of bins into fixed-length subsegments
    of ~target_length (default 0.005 m = 0.5 cm). All subsegments are
    target_length except the final one, which is shortened to end exactly
    at the last bin endpoint.

    Notes
    -----
    - Treats the input bins as a continuous polyline in (r, z).
    - If there is a discontinuity between bins (start != previous end),
      it closes the current chain and starts a new one (the "last segment"
      rule applies to each chain).
    - If ``fixed_vertices`` is given, a sub-bin boundary is forced at each
      of those (R, Z) points (the walking accumulator is reset there).
    - `solps_segments` is kept only for compatibility and is not used.
    - If your coordinates are in centimeters, call with `target_length=0.5`.

    Returns
    -------
    List[NewBin]
        New bins with lengths ~target_length, except the last of each chain.
    """
    # Default to 0.5 cm in meters
    if target_length is None:
        target_length = 0.001

    EPS = 1e-12

    if not bins:
        return []

    _fixed = list(fixed_vertices) if fixed_vertices else []
    _fixed_tol = 1e-6  # matching tolerance in metres

    def _is_fixed_vertex(r: float, z: float) -> bool:
        for fr, fz in _fixed:
            if abs(r - fr) < _fixed_tol and abs(z - fz) < _fixed_tol:
                return True
        return False

    def _emit(new_bins_list, bid, r0, z0, r1, z1):
        seglen = seg_length(r0, z0, r1, z1)
        if seglen > EPS:
            new_bins_list.append(NewBin(
                bin_id=bid,
                r_start=r0, z_start=z0,
                r_end=r1,   z_end=z1,
                arc_length=seglen,
            ))
            return bid + 1
        return bid

    new_bins: List[NewBin] = []
    bid = 1

    # Initialize the chain at the very first point
    chain_start_r = bins[0].r_start
    chain_start_z = bins[0].z_start
    last_boundary_r = chain_start_r  # start point of current output piece
    last_boundary_z = chain_start_z
    cursor_r = chain_start_r         # current walking point along the polyline
    cursor_z = chain_start_z
    accum_len = 0.0                  # length since last emitted boundary

    for b in bins:
        seg_start_r, seg_start_z = b.r_start, b.z_start
        seg_end_r, seg_end_z = b.r_end, b.z_end

        # If this bin starts away from our cursor, close leftover piece and
        # start a new chain at this bin's start.
        if seg_length(cursor_r, cursor_z, seg_start_r, seg_start_z) > EPS:
            # Emit leftover (short) segment of the previous chain, if any
            if seg_length(last_boundary_r, last_boundary_z, cursor_r, cursor_z) > EPS:
                bid = _emit(new_bins, bid, last_boundary_r, last_boundary_z, cursor_r, cursor_z)
            # Reset the chain at the start of the new bin
            last_boundary_r, last_boundary_z = seg_start_r, seg_start_z
            cursor_r, cursor_z = seg_start_r, seg_start_z
            accum_len = 0.0
        elif _is_fixed_vertex(seg_start_r, seg_start_z):
            # Force a sub-bin boundary at a fixed vertex so it is preserved
            # through merging later.
            if seg_length(last_boundary_r, last_boundary_z, cursor_r, cursor_z) > EPS:
                bid = _emit(new_bins, bid, last_boundary_r, last_boundary_z, cursor_r, cursor_z)
            last_boundary_r, last_boundary_z = seg_start_r, seg_start_z
            cursor_r, cursor_z = seg_start_r, seg_start_z
            accum_len = 0.0

        # Walk along this bin, emitting fixed-length segments as needed
        while True:
            # Remaining distance from cursor to end of this bin
            rem = seg_length(cursor_r, cursor_z, seg_end_r, seg_end_z)
            if rem < EPS:
                # We're effectively at the end of this bin
                break

            need = target_length - accum_len
            # Protect against numeric drift (e.g., accum_len ~ target_length)
            if need <= EPS:
                need = target_length

            if rem + accum_len < target_length - EPS:
                # Not enough length in this bin to reach the next boundary
                cursor_r, cursor_z = seg_end_r, seg_end_z
                accum_len += rem
                break

            # We can reach (or cross) the next boundary within this bin:
            # advance by 'need' along the current remainder
            frac = need / rem  # 0 < frac <= 1
            new_r = cursor_r + frac * (seg_end_r - cursor_r)
            new_z = cursor_z + frac * (seg_end_z - cursor_z)

            # Emit piece from last boundary to the new boundary
            bid = _emit(new_bins, bid, last_boundary_r, last_boundary_z, new_r, new_z)

            # Update state: new boundary becomes the start of next piece
            last_boundary_r, last_boundary_z = new_r, new_z
            cursor_r, cursor_z = new_r, new_z
            accum_len = 0.0  # reset since we've just hit a boundary

        # Proceed to next bin (cursor is at its end or at a boundary inside it)

    # End of all bins: if there's leftover distance since last boundary, emit it
    if seg_length(last_boundary_r, last_boundary_z, cursor_r, cursor_z) > EPS:
        bid = _emit(new_bins, bid, last_boundary_r, last_boundary_z, cursor_r, cursor_z)

    return new_bins


# ──────────────────────── Step 5: Map SOLPS data to bins ─────────────────────

import math
import numpy as np
from typing import List, Optional, Sequence, Tuple

# --- Attribute helpers (tolerant to different field names) --------------------

def _get_attr(obj, candidates, required=True, obj_name="object"):
    for name in candidates:
        if hasattr(obj, name):
            return getattr(obj, name)
    if required:
        raise AttributeError(f"{obj_name} is missing required attribute among: {candidates}")
    return None

def _get_bin_coords(bn):
    r0 = float(_get_attr(bn, ["r_start", "r0", "R_start", "R0", "rA"], obj_name="bin"))
    z0 = float(_get_attr(bn, ["z_start", "z0", "Z_start", "Z0", "zA"], obj_name="bin"))
    r1 = float(_get_attr(bn, ["r_end", "r1", "R_end", "R1", "rB"], obj_name="bin"))
    z1 = float(_get_attr(bn, ["z_end", "z1", "Z_end", "Z1", "zB"], obj_name="bin"))
    return r0, z0, r1, z1

def _get_seg_coords(seg):
    r0 = float(_get_attr(seg, ["r_start", "r0", "R_start", "R0", "rA"], obj_name="segment"))
    z0 = float(_get_attr(seg, ["z_start", "z0", "Z_start", "Z0", "zA"], obj_name="segment"))
    r1 = float(_get_attr(seg, ["r_end", "r1", "R_end", "R1", "rB"], obj_name="segment"))
    z1 = float(_get_attr(seg, ["z_end", "z1", "Z_end", "Z1", "zB"], obj_name="segment"))
    return r0, z0, r1, z1

def _get_seg_data(seg):
    data = getattr(seg, "data", None)
    return data if isinstance(data, dict) else {}

def _infer_columns_from_segments(segments) -> Sequence[str]:
    # Union of keys across segments to be safe if you don't pass columns
    seen, cols = set(), []
    for s in segments:
        for k in _get_seg_data(s).keys():
            if k not in seen:
                seen.add(k); cols.append(k)
    return cols

# --- Geometry helper: point-to-segments squared distances ---------------------

def _point_to_segments_dist2(P: np.ndarray, S0: np.ndarray, S1: np.ndarray, eps: float = 1e-15) -> np.ndarray:
    """
    Squared distance from point P[2] to each segment S0[i]->S1[i].
    S0, S1: (N,2). Returns (N,) squared distances.
    Handles degenerate segments safely.
    """
    v = S1 - S0                  # (N,2)
    w = P[None, :] - S0          # (N,2)
    vv = (v[:, 0]**2 + v[:, 1]**2)  # (N,)
    deg = vv <= eps              # zero-length segments
    t = np.zeros_like(vv)
    t[~deg] = (w[~deg, 0]*v[~deg, 0] + w[~deg, 1]*v[~deg, 1]) / (vv[~deg] + eps)
    t = np.clip(t, 0.0, 1.0)
    proj = S0 + t[:, None] * v
    d = P[None, :] - proj
    return d[:, 0]**2 + d[:, 1]**2

# --- Main mapping -------------------------------------------------------------

def map_solps_to_bins(
    bins: List["NewBin"],
    solps_segments: List["SOLPSSegment"],
    *,
    columns: Optional[Sequence[str]] = None,      # e.g. WLLD_COLUMNS; if None, try global; else infer
    max_dist: Optional[float] = None,             # cutoff in same units as R,Z; None => always assign
    fill_missing: float = float("nan"),           # fill if a column missing in a chosen segment
    n_subsample: Optional[int] = None,   # kept for compatibility (unused)
    return_debug: bool = False                    # return (idx, dist) arrays for inspection
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """
    For each bin:
      - Compute bin midpoint M_b = ((r_start + r_end)/2, (z_start + z_end)/2).
      - Compute the exact minimum distance from M_b to every SOLPS segment (continuous).
      - Choose the segment with minimal distance and copy its constant values.

    Notes
    -----
    - No interpolation: segments are assumed to carry one constant value per magnitude.
    - Deterministic tie-breaking via np.argmin (lowest index).
    - Bins/segments with invalid coords (NaN/Inf) are skipped or marked unassigned.
    - If `max_dist` is set and exceeded, the bin is left unassigned (NaN for requested columns).
    """
    if not bins:
        return None

    # Resolve columns
    if columns is not None:
        resolved_columns = list(columns)
    else:
        # Try global WLLD_COLUMNS if present; else infer from segments
        try:
            resolved_columns = list(WLLD_COLUMNS)  # type: ignore[name-defined]
        except NameError:
            resolved_columns = _infer_columns_from_segments(solps_segments)

    # Build arrays of valid SOLPS segments (filter out bad coords)
    S0_list, S1_list, valid_segments = [], [], []
    for seg in solps_segments:
        try:
            r0, z0, r1, z1 = _get_seg_coords(seg)
            if not (np.isfinite(r0) and np.isfinite(z0) and np.isfinite(r1) and np.isfinite(z1)):
                continue
            S0_list.append([r0, z0])
            S1_list.append([r1, z1])
            valid_segments.append(seg)
        except Exception:
            # Skip malformed segment
            continue

    if not valid_segments:
        # No usable segments; mark bins as unassigned
        for bn in bins:
            if resolved_columns:
                bn.data = {c: fill_missing for c in resolved_columns}
            bn.n_solps_contributors = 0
        return None

    S0 = np.asarray(S0_list, dtype=float)  # (N,2)
    S1 = np.asarray(S1_list, dtype=float)  # (N,2)

    chosen_indices = np.full(len(bins), -1, dtype=int)
    distances = np.full(len(bins), np.nan, dtype=float)

    for i, bn in enumerate(bins):
        try:
            r0, z0, r1, z1 = _get_bin_coords(bn)
            if not (np.isfinite(r0) and np.isfinite(z0) and np.isfinite(r1) and np.isfinite(z1)):
                if resolved_columns:
                    bn.data = {c: fill_missing for c in resolved_columns}
                bn.n_solps_contributors = 0
                continue

            # Bin midpoint
            Mb = np.array([(r0 + r1)/2.0, (z0 + z1)/2.0], dtype=float)

            # Distances from midpoint to each SOLPS segment (continuous)
            d2 = _point_to_segments_dist2(Mb, S0, S1)
            j = int(np.argmin(d2))
            dist = float(math.sqrt(d2[j]))

            # Optional cutoff
            if (max_dist is not None) and (dist > max_dist):
                if resolved_columns:
                    bn.data = {c: fill_missing for c in resolved_columns}
                bn.n_solps_contributors = 0
                chosen_indices[i] = -1
                distances[i] = dist
                continue

            # Assign chosen segment's constant values
            seg = valid_segments[j]
            seg_data = _get_seg_data(seg)
            if resolved_columns:
                bn.data = {c: float(seg_data.get(c, fill_missing)) for c in resolved_columns}
            else:
                # If no explicit columns, copy the whole segment dict
                bn.data = dict(seg_data)
            bn.n_solps_contributors = 1

            chosen_indices[i] = j
            distances[i] = dist

        except Exception:
            # On any unexpected failure, mark bin unassigned
            if resolved_columns:
                bn.data = {c: fill_missing for c in resolved_columns}
            bn.n_solps_contributors = 0
            chosen_indices[i] = -1
            distances[i] = np.nan

    return (chosen_indices, distances) if return_debug else None

# ──────────────────────── Step 6b: Merge similar adjacent bins ───────────────

def merge_similar_bins(
    bins: List[NewBin],
    *,
    rtol: float = 0.15,
    rtol_hot: float = 0.08,
    hot_thresh: float = 1e6,
    atol_wtot: float = 0.0,
    min_length: float = 0.0,
    max_length: float = 0.20,
    connect_tol: float = 1e-6,
    fixed_vertices: Optional[List[Tuple[float, float]]] = None,
) -> List[NewBin]:
    """
    Merge adjacent sub-bins using a greedy left-to-right sweep.

    Algorithm (greedy-extend):
      Start a new bin with the first sub-bin.  Keep adding the next
      sub-bin as long as:
        • the shared vertex is not a fixed vertex,
        • the bins are connected (within connect_tol),
        • the merged length does not exceed max_length,
        • the merged average Wtot stays within *rtol* (or *rtol_hot* if any
          constituent has Wtot > hot_thresh) of every constituent sub-bin,
        • the merged average D-flux stays within *rtol* of every constituent.
      When a sub-bin fails the check, finalize the current merged bin and
      start a new one.

    The tolerance check is performed at a coarser resolution (CHECK_RES)
    rather than at every 1 mm sub-bin, so that natural SOLPS gradients
    within a single cell do not create artificial bin boundaries.

    fixed_vertices : list of (R, Z) points that must remain as bin
        boundaries (e.g. sharp corners created by the snap/trim step).
    """
    if not bins:
        return []

    CHECK_RES = 0.08  # 80 mm – coarser check resolution for tolerance

    _fixed = []
    _fixed_tol = 1e-4
    if fixed_vertices:
        _fixed = list(fixed_vertices)

    def _is_fixed(r: float, z: float) -> bool:
        for fr, fz in _fixed:
            if abs(r - fr) < _fixed_tol and abs(z - fz) < _fixed_tol:
                return True
        return False

    def _dflux(b: NewBin) -> float:
        return b.data.get("flxi_D", 0.0) + b.data.get("flxa_D", 0.0)

    def _connected(a_end_r, a_end_z, b_start_r, b_start_z) -> bool:
        return seg_length(a_end_r, a_end_z, b_start_r, b_start_z) <= connect_tol

    def _make_merged_bin(grp: List[NewBin]) -> NewBin:
        """Create a single NewBin from a group of sub-bins."""
        total_len = sum(b.arc_length for b in grp)
        merged_data: Dict[str, float] = {}
        all_keys: set = set()
        for b in grp:
            all_keys |= set(b.data.keys())
        for k in all_keys:
            if total_len > 0:
                merged_data[k] = sum(b.data.get(k, 0.0) * b.arc_length for b in grp) / total_len
            else:
                merged_data[k] = sum(b.data.get(k, 0.0) for b in grp) / max(len(grp), 1)
        n_contrib = sum(b.n_solps_contributors for b in grp)
        return NewBin(
            bin_id=grp[0].bin_id,
            r_start=grp[0].r_start, z_start=grp[0].z_start,
            r_end=grp[-1].r_end,    z_end=grp[-1].z_end,
            arc_length=total_len,
            data=merged_data,
            n_solps_contributors=n_contrib,
        )

    # --- Pre-aggregate sub-bins into coarser "check-bins" (CHECK_RES) --------
    # These coarser bins are used for tolerance checks only; the final merged
    # bin values still use all constituent 1 mm sub-bins for accuracy.
    check_bins: List[NewBin] = []
    i = 0
    while i < len(bins):
        grp = [bins[i]]
        acc_len = bins[i].arc_length
        j = i + 1
        while j < len(bins) and acc_len < CHECK_RES:
            if not _connected(bins[j-1].r_end, bins[j-1].z_end,
                              bins[j].r_start, bins[j].z_start):
                break
            if _is_fixed(bins[j-1].r_end, bins[j-1].z_end):
                break
            grp.append(bins[j])
            acc_len += bins[j].arc_length
            j += 1
        check_bins.append(_make_merged_bin(grp))
        i = j

    # --- Greedy-extend merge -------------------------------------------------
    # Work with check_bins; for each candidate merged bin, verify that the
    # merged average stays within tolerance of every constituent check-bin.
    result_groups: List[List[NewBin]] = []  # each element = list of check_bins
    current: List[NewBin] = [check_bins[0]]

    for k in range(1, len(check_bins)):
        cb = check_bins[k]
        prev = current[-1]
        # Structural checks
        can_add = True
        if not _connected(prev.r_end, prev.z_end, cb.r_start, cb.z_start):
            can_add = False
        elif _is_fixed(prev.r_end, prev.z_end):
            can_add = False
        else:
            # Length check
            cand = current + [cb]
            total_len = sum(b.arc_length for b in cand)
            if total_len > max_length:
                can_add = False
            else:
                # Tolerance check: merged avg vs each constituent check-bin
                avg_w = sum(b.data.get("Wtot", 0.0) * b.arc_length for b in cand) / total_len
                avg_d = sum(_dflux(b) * b.arc_length for b in cand) / total_len

                is_hot = any(b.data.get("Wtot", 0.0) > hot_thresh for b in cand)
                w_tol = rtol_hot if is_hot else rtol

                for b in cand:
                    bw = b.data.get("Wtot", 0.0)
                    bd = _dflux(b)
                    if abs(avg_w - bw) / max(abs(avg_w), abs(bw), 1.0) > w_tol:
                        can_add = False
                        break
                    if abs(avg_d - bd) / max(abs(avg_d), abs(bd), 1.0) > rtol:
                        can_add = False
                        break

        if can_add:
            current.append(cb)
        else:
            result_groups.append(current)
            current = [cb]
    result_groups.append(current)

    # Build final merged bins – use the original 1 mm sub-bins for accuracy.
    # Map each check-bin back to its arc-length span in the original bins list.
    # Actually, each result_group is a list of check_bins (already merged from
    # 1 mm sub-bins).  We can merge the check_bins in each group.
    out = [_make_merged_bin(grp) for grp in result_groups]

    # --- Absorb tiny edge bins into their neighbor -------------------------
    min_edge = MIN_EDGE_BIN_LENGTH
    changed = True
    while changed and len(out) >= 2:
        changed = False
        # First bin too short → merge with second
        if out[0].arc_length < min_edge and _connected(
                out[0].r_end, out[0].z_end, out[1].r_start, out[1].z_start):
            out[0:2] = [_make_merged_bin(out[0:2])]
            changed = True
        # Last bin too short → merge with second-to-last
        if len(out) >= 2 and out[-1].arc_length < min_edge and _connected(
                out[-2].r_end, out[-2].z_end, out[-1].r_start, out[-1].z_start):
            out[-2:] = [_make_merged_bin(out[-2:])]
            changed = True

    # Renumber
    for k, b in enumerate(out, start=1):
        b.bin_id = k

    return out


# ──────────────────────── FW Bin Mapping ─────────────────────────────────────

def parse_fw_bins_from_csv(path: str) -> List[NewBin]:
    """Read First-Wall bins from input_table.csv (location == 'FW').

    Returns one NewBin per unique bin number, ordered by bin number.
    Coordinates are taken from the first row for each bin.
    """
    bins_dict: Dict[int, NewBin] = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            loc = row["location"].strip()
            if loc != "FW":
                continue
            bn = int(row["Bin number"].strip())
            if bn in bins_dict:
                continue
            z_start = float(row["Z_start (m)"].strip())
            r_start = float(row["R_start (m)"].strip())
            z_end   = float(row["Z_end (m)"].strip())
            r_end   = float(row["R_end (m)"].strip())
            arc = seg_length(r_start, z_start, r_end, z_end)
            bins_dict[bn] = NewBin(
                bin_id=bn, r_start=r_start, z_start=z_start,
                r_end=r_end, z_end=z_end, arc_length=arc,
            )
    return [bins_dict[k] for k in sorted(bins_dict)]


def parse_fw_bins_from_coords(path: str) -> List[NewBin]:
    """Read First-Wall bins from fw_bins_coords.dat.

    File format (whitespace-separated):
        # bin_number Z_start(m) R_start(m) Z_end(m) R_end(m) arc_length(m)
        1  -2.51000  4.10000  -1.50000  4.10000  1.01000
        ...

    Returns list of NewBin objects ordered by bin number.
    """
    bins: List[NewBin] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                bid = int(parts[0])
                z_start = float(parts[1])
                r_start = float(parts[2])
                z_end   = float(parts[3])
                r_end   = float(parts[4])
                arc     = float(parts[5])
            except (ValueError, IndexError):
                continue
            bins.append(NewBin(
                bin_id=bid, r_start=r_start, z_start=z_start,
                r_end=r_end, z_end=z_end, arc_length=arc,
            ))
    return sorted(bins, key=lambda b: b.bin_id)


def map_solps_to_fw_bins(
    fw_bins: List[NewBin],
    solps_segments: List[SOLPSSegment],
    subseg_len: float = 0.001,   # 1 mm sub-segments
) -> None:
    """Map SOLPS fluxes onto FW bins via length-weighted averaging.

    Strategy:
      1. Build array of FW bin edges (segment start→end).
      2. For every SOLPS segment, chop it into ~1 mm sub-segments.
      3. For each sub-segment midpoint, find the closest FW bin (point-to-segment).
      4. Accumulate (value × sub-segment length) and total length per FW bin.
      5. Final value = sum(value × length) / sum(length) for each column.
    """
    import math as _math

    n_fw = len(fw_bins)
    if n_fw == 0:
        return

    # FW bin edges as arrays for vectorised distance computation
    FW_S0 = np.array([[b.r_start, b.z_start] for b in fw_bins], dtype=float)  # (n_fw, 2)
    FW_S1 = np.array([[b.r_end,   b.z_end]   for b in fw_bins], dtype=float)  # (n_fw, 2)

    # Accumulators
    columns = list(WLLD_COLUMNS)
    weighted_sums = np.zeros((n_fw, len(columns)), dtype=float)
    total_length  = np.zeros(n_fw, dtype=float)
    n_contrib     = np.zeros(n_fw, dtype=int)

    for seg in solps_segments:
        # Only process segments that have data
        if not seg.data:
            continue
        vals = np.array([seg.data.get(c, 0.0) for c in columns], dtype=float)

        # Segment geometry
        dr = seg.r_end - seg.r_start
        dz = seg.z_end - seg.z_start
        L  = _math.sqrt(dr*dr + dz*dz)
        if L < 1e-12:
            continue

        # Number of sub-segments
        n_sub = max(1, int(round(L / subseg_len)))
        sub_L = L / n_sub

        for k in range(n_sub):
            t = (k + 0.5) / n_sub
            mid = np.array([seg.r_start + t * dr,
                            seg.z_start + t * dz], dtype=float)

            # Distance from midpoint to each FW bin edge
            d2 = _point_to_segments_dist2(mid, FW_S0, FW_S1)
            j = int(np.argmin(d2))

            weighted_sums[j] += vals * sub_L
            total_length[j]  += sub_L
            n_contrib[j]     += 1

    # Compute weighted averages and store in bins
    for i, b in enumerate(fw_bins):
        if total_length[i] > 0:
            avg = weighted_sums[i] / total_length[i]
            b.data = {c: float(avg[k]) for k, c in enumerate(columns)}
        else:
            b.data = {c: 0.0 for c in columns}
        b.n_solps_contributors = int(n_contrib[i])


# ──────────────────────── Step 7: Output ─────────────────────────────────────

def write_bin_coords(bins: List[NewBin], path: str) -> None:
    """Write bin coordinates in the same format as divertor_bins_coords.dat."""
    with open(path, "w") as fh:
        fh.write("# bin_number Z_start(m) R_start(m) Z_end(m) R_end(m) "
                 "arc_length(m)\n")
        for bn in bins:
            fh.write(f"{bn.bin_id:4d}  {bn.z_start:10.5f}  {bn.r_start:10.5f}"
                     f"  {bn.z_end:10.5f}  {bn.r_end:10.5f}"
                     f"  {bn.arc_length:10.5f}\n")
    print(f"Wrote {len(bins)} bin coordinates to {path}")


def write_bin_data(bins: List[NewBin], path: str) -> None:
    """Write full mapped SOLPS data for each bin."""
    with open(path, "w") as fh:
        header = "# bin_number  " + "  ".join(WLLD_COLUMNS) + "  n_contributors\n"
        fh.write(header)
        for bn in bins:
            vals = "  ".join(f"{bn.data.get(c, 0.0):14.6E}" for c in WLLD_COLUMNS)
            fh.write(f"{bn.bin_id:4d}  {vals}  {bn.n_solps_contributors}\n")
    print(f"Wrote bin data to {path}")


def write_binned_flux_csv(div_bins: List[NewBin],
                          fw_bins: List[NewBin],
                          path: str) -> None:
    """Write merged FW + DIV bins in Binned_Flux_Data.dat CSV format.

    Column order:
        Bin_Index, Flux_Ion, Flux_Atom, E_ion, E_atom,
        alpha_ion, alpha_atom, heat_total, heat_ion

    FW bins are numbered first (starting at 0), then DIV bins continue
    the index sequentially.
    """
    ALPHA_ION = 60.0
    ALPHA_ATOM = 45.0
    header = ("Bin_Index,Flux_Ion,Flux_Atom,E_ion,E_atom,"
              "alpha_ion,alpha_atom,heat_total,heat_ion\n")

    with open(path, "w") as fh:
        fh.write(header)
        idx = 0
        # FW bins first
        for b in sorted(fw_bins, key=lambda x: x.bin_id):
            fh.write(f"{idx},"
                     f"{b.data.get('flxi_D', 0.0):.6E},"
                     f"{b.data.get('flxa_D', 0.0):.6E},"
                     f"{b.data.get('Eavi_D', 0.0):.6E},"
                     f"{b.data.get('Eava_D', 0.0):.6E},"
                     f"{ALPHA_ION:.1f},"
                     f"{ALPHA_ATOM:.1f},"
                     f"{b.data.get('Wtot', 0.0):.6E},"
                     f"{b.data.get('Wpls', 0.0):.6E}\n")
            idx += 1
        # DIV bins next
        for b in sorted(div_bins, key=lambda x: x.bin_id):
            fh.write(f"{idx},"
                     f"{b.data.get('flxi_D', 0.0):.6E},"
                     f"{b.data.get('flxa_D', 0.0):.6E},"
                     f"{b.data.get('Eavi_D', 0.0):.6E},"
                     f"{b.data.get('Eava_D', 0.0):.6E},"
                     f"{ALPHA_ION:.1f},"
                     f"{ALPHA_ATOM:.1f},"
                     f"{b.data.get('Wtot', 0.0):.6E},"
                     f"{b.data.get('Wpls', 0.0):.6E}\n")
            idx += 1
    print(f"Wrote {idx} bins ({len(fw_bins)} FW + {len(div_bins)} DIV) "
          f"to {path}")


# ──────────────────────── Step 8: Profile Plots ──────────────────────────────

def _compute_arc_positions(bins: List[NewBin]) -> np.ndarray:
    """Return array of cumulative arc-length midpoints for a list of bins.

    Bins are assumed to be ordered; the returned positions represent the
    midpoint of each bin along the arc.
    """
    positions = np.zeros(len(bins))
    cumulative = 0.0
    for i, b in enumerate(bins):
        positions[i] = cumulative + b.arc_length / 2.0
        cumulative += b.arc_length
    return positions


def _project_solps_onto_bins(
    solps_segs: List[SOLPSSegment],
    bins_sorted: List[NewBin],
    max_dist: float = 0.15,
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """Project SOLPS segment midpoints onto the bin arc-length axis.

    For each SOLPS segment whose midpoint is within *max_dist* of a bin
    midpoint, assign it the arc-length position of that nearest bin.

    Returns (positions, data_dict)  –  both sorted by arc-length.
    *data_dict* maps WLLD column names to value arrays.
    """
    if not bins_sorted or not solps_segs:
        return np.array([]), {}

    n_bins = len(bins_sorted)
    bin_mids = np.array(
        [[0.5 * (b.r_start + b.r_end), 0.5 * (b.z_start + b.z_end)]
         for b in bins_sorted], dtype=float,
    )  # (n_bins, 2)

    cum_arc = np.zeros(n_bins)
    for i in range(1, n_bins):
        cum_arc[i] = cum_arc[i - 1] + bins_sorted[i - 1].arc_length

    positions: List[float] = []
    col_lists: Dict[str, List[float]] = {c: [] for c in WLLD_COLUMNS}

    for seg in solps_segs:
        mid = np.array([seg.r_centre, seg.z_centre])
        dists = np.linalg.norm(bin_mids - mid, axis=1)
        best = int(np.argmin(dists))
        if dists[best] > max_dist:
            continue
        pos = cum_arc[best] + bins_sorted[best].arc_length / 2.0
        positions.append(pos)
        for c in WLLD_COLUMNS:
            col_lists[c].append(seg.data.get(c, 0.0))

    if not positions:
        return np.array([]), {}

    order = np.argsort(positions)
    pos_arr = np.array(positions)[order]
    data_dict = {c: np.array(col_lists[c])[order] for c in WLLD_COLUMNS}
    return pos_arr, data_dict


def plot_profiles(div_bins: List[NewBin],
                  fw_bins: List[NewBin],
                  out_dir: str,
                  solps_all: Optional[List[SOLPSSegment]] = None) -> None:
    """Create arc-length profile plots for fluxes, energies, and heat loads.

    Produces three PNG files in *out_dir*:
        profile_flux.png   – Ion and atom particle flux along arc-length
        profile_energy.png – Ion and atom mean energy along arc-length
        profile_heat.png   – Total heat load and ion heat along arc-length

    Each plot has a 2×2 grid:
        top-left: FW (left quantity)    top-right: FW (right quantity)
        bot-left: DIV (left quantity)   bot-right: DIV (right quantity)

    Y-axis is always log-scale.  Raw SOLPS values are shown as a light-grey
    background line.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    div_sorted = sorted(div_bins, key=lambda b: b.bin_id)
    fw_sorted  = sorted(fw_bins,  key=lambda b: b.bin_id)

    div_arc = _compute_arc_positions(div_sorted)
    fw_arc  = _compute_arc_positions(fw_sorted)

    # Project raw SOLPS segments onto each region's arc-length axis
    solps_div_arc, solps_div_data = _project_solps_onto_bins(
        solps_all or [], div_sorted)
    solps_fw_arc, solps_fw_data = _project_solps_onto_bins(
        solps_all or [], fw_sorted)

    def _get(bins, key):
        return np.array([b.data.get(key, 0.0) for b in bins])

    def _make_2x2(title_left, title_right, key_left, key_right,
                  ylabel_left, ylabel_right, fname):
        fig, axes = plt.subplots(2, 2, figsize=(16, 10))
        fig.suptitle("Profiles along arc-length", fontsize=14, y=0.98)

        datasets = [
            (fw_sorted, fw_arc, solps_fw_arc, solps_fw_data, "FW", 0),
            (div_sorted, div_arc, solps_div_arc, solps_div_data, "DIV", 1),
        ]
        for bins_list, arc, s_arc, s_data, region, row in datasets:
            for col_idx, (key, title, ylabel, color) in enumerate([
                (key_left, title_left, ylabel_left, "tab:blue"),
                (key_right, title_right, ylabel_right, "tab:red"),
            ]):
                ax = axes[row, col_idx]

                # Raw SOLPS background (light grey)
                if s_data and key in s_data and len(s_arc) > 0:
                    s_vals = np.where(s_data[key] > 0, s_data[key], np.nan)
                    ax.plot(s_arc, s_vals, color="lightgrey",
                            linewidth=0.8, zorder=1, label="SOLPS raw")

                # Binned values (step plot)
                vals = _get(bins_list, key)
                vals_safe = np.where(vals > 0, vals, np.nan)
                ax.step(arc, vals_safe, where="mid", linewidth=1.5,
                        color=color, zorder=2, label="Binned")

                ax.set_yscale("log")
                ax.yaxis.set_major_formatter(
                    mticker.FuncFormatter(
                        lambda x, _: f"{x:.1E}" if x > 0 else "0"))
                ax.set_title(f"{region} – {title}", fontsize=11)
                ax.set_xlabel("Arc-length (m)")
                ax.set_ylabel(ylabel)
                ax.grid(True, linewidth=0.3, alpha=0.5)
                ax.legend(fontsize=8)

        plt.tight_layout()
        out_path = os.path.join(out_dir, fname)
        plt.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved profile plot: {out_path}")

    # 1. Flux profiles: ion flux | atom flux
    _make_2x2("Ion flux (flxi_D)", "Atom flux (flxa_D)",
              "flxi_D", "flxa_D",
              "Particle flux (1/m²/s)", "Particle flux (1/m²/s)",
              "profile_flux.png")

    # 2. Energy profiles: ion energy | atom energy
    _make_2x2("Ion energy (Eavi_D)", "Atom energy (Eava_D)",
              "Eavi_D", "Eava_D",
              "Mean energy (eV)", "Mean energy (eV)",
              "profile_energy.png")

    # 3. Heat load profiles: total heat | ion heat
    _make_2x2("Total heat load (Wtot)", "Ion heat (Wpls)",
              "Wtot", "Wpls",
              "Heat load (W/m²)", "Heat load (W/m²)",
              "profile_heat.png")


# ──────────────────────── Step 9: Diagnostic Plotting ────────────────────────

def plot_all(solps_all: List[SOLPSSegment],
             solps_filtered: List[SOLPSSegment],
             div_groups: List[List[Tuple[float, float]]],
             bins: List[NewBin],
             path: str,
             conservation: Optional[Tuple[float, float]] = None,
             fw_bins: Optional[List[NewBin]] = None) -> None:
    """Create a combined diagnostic plot.

    conservation : (pct_wtot, pct_dflux) – if given, printed on the figure.
    fw_bins      : First-Wall bins with mapped SOLPS data.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    import matplotlib.cm as cm

    has_fw = fw_bins is not None and len(fw_bins) > 0
    nrows = 2 if has_fw else 1
    fig, axes = plt.subplots(nrows, 3, figsize=(24, 10 * nrows))
    if nrows == 1:
        axes = axes[np.newaxis, :]  # make 2-D so indexing works

    # Compute divertor bounding box from bins (with margin) for all panels
    if bins:
        _br = [b.r_start for b in bins] + [b.r_end for b in bins]
        _bz = [b.z_start for b in bins] + [b.z_end for b in bins]
        _margin = 0.15
        _r_lo, _r_hi = min(_br) - _margin, max(_br) + _margin
        _z_lo, _z_hi = min(_bz) - _margin, max(_bz) + _margin
    else:
        _r_lo, _r_hi, _z_lo, _z_hi = 3.5, 7.0, -5.0, -2.0

    # ── Row 0, Panel 1: Overview (divertor region only) ──
    ax = axes[0, 0]
    ax.set_title("Divertor – all layers superposed", fontsize=11)

    # SOLPS segments (all in divertor box, grey, thin)
    for seg in solps_all:
        if _r_lo <= seg.r_centre <= _r_hi and _z_lo <= seg.z_centre <= _z_hi:
            ax.plot([seg.r_start, seg.r_end], [seg.z_start, seg.z_end],
                    color="silver", linewidth=0.4, zorder=1)

    # SOLPS segments (filtered, blue)
    for seg in solps_filtered:
        ax.plot([seg.r_start, seg.r_end], [seg.z_start, seg.z_end],
                color="royalblue", linewidth=0.8, zorder=2)

    # Divertor polylines (red)
    for grp in div_groups:
        rs = [p[0] for p in grp]
        zs = [p[1] for p in grp]
        ax.plot(rs, zs, color="red", linewidth=1.8, marker=".", markersize=3,
                zorder=3, label="Real divertor" if grp is div_groups[0] else "")

    # New bins (green, thick)
    for i, bn in enumerate(bins):
        ax.plot([bn.r_start, bn.r_end], [bn.z_start, bn.z_end],
                color="limegreen", linewidth=2.5, solid_capstyle="round",
                zorder=4, label="New bins" if i == 0 else "")
        # Mark endpoints
        ax.plot(bn.r_start, bn.z_start, "g.", markersize=3, zorder=5)
        ax.plot(bn.r_end, bn.z_end, "g.", markersize=3, zorder=5)

    ax.set_xlabel("R (m)")
    ax.set_ylabel("Z (m)")
    ax.set_aspect("equal", "box")
    ax.set_xlim(_r_lo, _r_hi)
    ax.set_ylim(_z_lo, _z_hi)
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, linewidth=0.3, alpha=0.5)

    # ── Row 0, Panel 2: Bins coloured by total heat flux (log scale) ──
    ax = axes[0, 1]
    ax.set_title("New bins – Wtot (W/m²)  [log]", fontsize=11)
    if bins:
        from matplotlib.colors import LogNorm
        wtot_vals = np.array([bn.data.get("Wtot", 0.0) for bn in bins])
        # Floor zeros/negatives to a small positive value for log scale
        wtot_safe = np.where(wtot_vals > 0, wtot_vals, np.nan)
        vmin_w = np.nanmin(wtot_safe) if np.any(~np.isnan(wtot_safe)) else 1.0
        vmax_w = np.nanmax(wtot_safe) if np.any(~np.isnan(wtot_safe)) else 1.0
        # Replace non-positive with vmin for display
        wtot_plot = np.where(wtot_vals > 0, wtot_vals, vmin_w)
        norm = LogNorm(vmin=vmin_w, vmax=vmax_w)
        cmap = cm.plasma

        segs_lc = [[(bn.r_start, bn.z_start), (bn.r_end, bn.z_end)]
                    for bn in bins]
        lc = LineCollection(segs_lc, cmap=cmap, norm=norm, linewidths=3)
        lc.set_array(wtot_plot)
        ax.add_collection(lc)
        fig.colorbar(lc, ax=ax, label="Wtot (W/m²)", shrink=0.7)

    # Divertor outline
    for grp in div_groups:
        rs = [p[0] for p in grp]
        zs = [p[1] for p in grp]
        ax.plot(rs, zs, "r--", linewidth=0.8, alpha=0.5)

    ax.set_xlabel("R (m)")
    ax.set_ylabel("Z (m)")
    ax.set_aspect("equal", "box")
    ax.grid(True, linewidth=0.3, alpha=0.5)
    ax.autoscale_view()

    # ── Row 0, Panel 3: Bins coloured by total D flux = flxi_D + flxa_D (log scale) ──
    ax = axes[0, 2]
    ax.set_title("New bins – D flux (ion+atom) (1/m²/s)  [log]", fontsize=11)
    if bins:
        from matplotlib.colors import LogNorm
        flux_vals = np.array([
            bn.data.get("flxi_D", 0.0) + bn.data.get("flxa_D", 0.0)
            for bn in bins
        ])
        flux_safe = np.where(flux_vals > 0, flux_vals, np.nan)
        vmin_f = np.nanmin(flux_safe) if np.any(~np.isnan(flux_safe)) else 1.0
        vmax_f = np.nanmax(flux_safe) if np.any(~np.isnan(flux_safe)) else 1.0
        flux_plot = np.where(flux_vals > 0, flux_vals, vmin_f)
        norm = LogNorm(vmin=vmin_f, vmax=vmax_f)
        cmap = cm.inferno

        segs_lc = [[(bn.r_start, bn.z_start), (bn.r_end, bn.z_end)]
                    for bn in bins]
        lc = LineCollection(segs_lc, cmap=cmap, norm=norm, linewidths=3)
        lc.set_array(flux_plot)
        ax.add_collection(lc)
        fig.colorbar(lc, ax=ax, label="D flux (flxi_D + flxa_D) (1/m²/s)", shrink=0.7)

    for grp in div_groups:
        rs = [p[0] for p in grp]
        zs = [p[1] for p in grp]
        ax.plot(rs, zs, "r--", linewidth=0.8, alpha=0.5)

    ax.set_xlabel("R (m)")
    ax.set_ylabel("Z (m)")
    ax.set_aspect("equal", "box")
    ax.grid(True, linewidth=0.3, alpha=0.5)
    ax.autoscale_view()

    # ── Row 1: FW bins (full tokamak cross-section) ──
    if has_fw:
        from matplotlib.colors import LogNorm as _LN

        # Combine divertor bins + FW bins for shared colour scales
        all_bins = list(bins) + list(fw_bins)

        # Helper: collect all positive values for a field across bins
        def _pos(field_fn):
            v = np.array([field_fn(b) for b in all_bins])
            return v[v > 0]

        # -- Panel (1,0): overview of FW + divertor bins geometry --
        ax = axes[1, 0]
        ax.set_title("Full wall – FW & Divertor bins", fontsize=11)
        # SOLPS segments (all, grey)
        for seg in solps_all:
            ax.plot([seg.r_start, seg.r_end], [seg.z_start, seg.z_end],
                    color="silver", linewidth=0.3, zorder=1)
        # Divertor bins (green)
        for i, bn in enumerate(bins):
            ax.plot([bn.r_start, bn.r_end], [bn.z_start, bn.z_end],
                    color="limegreen", linewidth=2.0, solid_capstyle="round",
                    zorder=3, label="Div bins" if i == 0 else "")
        # FW bins (cyan)
        for i, bn in enumerate(fw_bins):
            ax.plot([bn.r_start, bn.r_end], [bn.z_start, bn.z_end],
                    color="cyan", linewidth=2.0, solid_capstyle="round",
                    zorder=3, label="FW bins" if i == 0 else "")
        ax.set_xlabel("R (m)"); ax.set_ylabel("Z (m)")
        ax.set_aspect("equal", "box")
        ax.legend(fontsize=8, loc="best")
        ax.grid(True, linewidth=0.3, alpha=0.5)
        ax.autoscale_view()

        # -- Panel (1,1): Wtot heatmap for FW + Div --
        ax = axes[1, 1]
        ax.set_title("All bins – Wtot (W/m²)  [log]", fontsize=11)
        pw = _pos(lambda b: b.data.get("Wtot", 0.0))
        if pw.size > 0:
            norm_w = _LN(vmin=pw.min(), vmax=pw.max())
            wtot_all = np.array([b.data.get("Wtot", 0.0) for b in all_bins])
            wtot_all = np.where(wtot_all > 0, wtot_all, pw.min())
            segs_lc = [[(b.r_start, b.z_start), (b.r_end, b.z_end)] for b in all_bins]
            lc = LineCollection(segs_lc, cmap=cm.plasma, norm=norm_w, linewidths=2.5)
            lc.set_array(wtot_all)
            ax.add_collection(lc)
            fig.colorbar(lc, ax=ax, label="Wtot (W/m²)", shrink=0.7)
        ax.set_xlabel("R (m)"); ax.set_ylabel("Z (m)")
        ax.set_aspect("equal", "box")
        ax.grid(True, linewidth=0.3, alpha=0.5)
        ax.autoscale_view()

        # -- Panel (1,2): D flux heatmap for FW + Div --
        ax = axes[1, 2]
        ax.set_title("All bins – D flux (ion+atom) (1/m²/s)  [log]", fontsize=11)
        pf = _pos(lambda b: b.data.get("flxi_D", 0.0) + b.data.get("flxa_D", 0.0))
        if pf.size > 0:
            norm_f = _LN(vmin=pf.min(), vmax=pf.max())
            flux_all = np.array([
                b.data.get("flxi_D", 0.0) + b.data.get("flxa_D", 0.0)
                for b in all_bins
            ])
            flux_all = np.where(flux_all > 0, flux_all, pf.min())
            segs_lc = [[(b.r_start, b.z_start), (b.r_end, b.z_end)] for b in all_bins]
            lc = LineCollection(segs_lc, cmap=cm.inferno, norm=norm_f, linewidths=2.5)
            lc.set_array(flux_all)
            ax.add_collection(lc)
            fig.colorbar(lc, ax=ax, label="D flux (1/m²/s)", shrink=0.7)
        ax.set_xlabel("R (m)"); ax.set_ylabel("Z (m)")
        ax.set_aspect("equal", "box")
        ax.grid(True, linewidth=0.3, alpha=0.5)
        ax.autoscale_view()

    plt.tight_layout()

    # Add conservation text if available
    if conservation is not None:
        pct_w, pct_d = conservation
        txt = (f"Conservation (bins / SOLPS):\n"
               f"  Wtot : {pct_w:.1f}%\n"
               f"  D flux: {pct_d:.1f}%")
        fig.text(0.01, 0.01, txt, fontsize=9, family="monospace",
                 va="bottom", ha="left",
                 bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="gray", alpha=0.85))

    plt.savefig(path, dpi=200, bbox_inches="tight")
    print(f"Saved combined plot to {path}")
    plt.close(fig)


# ──────────────────────── Main Pipeline ──────────────────────────────────────

def main() -> None:
    print("=" * 70)
    print("  Automatic Binning")
    print("=" * 70)

    # ── Validate input files exist ──
    print(f"\n  SOLPS data dir : {SOLPS_DIR}")
    print(f"  Fixed bins dir : {FIXED_BINS_DIR}")
    print(f"  Coords-to-bin  : {COORDS_TO_BIN_DIR}")
    print(f"  Output dir     : {OUTPUT_DIR}")
    for p, label in [(WLLD_PATH, "wlld"), (WLLY_PATH, "wlly")]:
        if not os.path.exists(p):
            sys.exit(f"ERROR: {label} file not found: {p}")

    # 1. Parse SOLPS input files
    print("\n[1/7] Parsing SOLPS input files ...")
    solps_segments = parse_wlly(WLLY_PATH)
    solps_data_by_x = parse_wlld(WLLD_PATH)
    print(f"  wlly segments  : {len(solps_segments)}")
    print(f"  wlld entries   : {len(solps_data_by_x)}")

    # 2. Link SOLPS geometry <-> data (match by x-coordinate)
    print("\n[2/7] Linking SOLPS geometry with plasma data (by x-value) ...")
    link_solps_data(solps_segments, solps_data_by_x)
    n_linked = sum(1 for s in solps_segments if s.data.get("Wtot", 0.0) != 0.0)
    print(f"  Segments with non-zero Wtot: {n_linked}/{len(solps_segments)}")

    # ══════════════════════════════════════════════════════════════════════
    # FIXED BINS  (Reactor_coordinates/fixed_bins/)
    # Treatment: parse pre-defined bins, map SOLPS values (no splitting/
    # shadow-filtering/merging).
    # ══════════════════════════════════════════════════════════════════════
    all_fixed_bins: List[NewBin] = []
    if FIXED_BIN_FILES:
        print(f"\n── Fixed bins ({len(FIXED_BIN_FILES)} file(s)) ──")
        for fb_path in FIXED_BIN_FILES:
            tag = os.path.splitext(os.path.basename(fb_path))[0]
            # Clean up tag: strip trailing '_coords' to avoid names like mapped_x_coords_coords
            if tag.endswith("_coords"):
                tag = tag[:-len("_coords")]
            print(f"\n  [{tag}] Parsing {os.path.basename(fb_path)} ...")
            fb = parse_fw_bins_from_coords(fb_path)
            print(f"  [{tag}] Bins read: {len(fb)}  "
                  f"(bins {fb[0].bin_id}–{fb[-1].bin_id})" if fb else
                  f"  [{tag}] WARNING: no bins found in file!")

            if fb:
                print(f"  [{tag}] Mapping SOLPS data (length-weighted average) ...")
                # Exclude divertor-region SOLPS segments to avoid contamination
                solps_for_fixed = [s for s in solps_segments
                                   if not (0.5*(s.z_start + s.z_end) < -2.5
                                           and 0.5*(s.r_start + s.r_end) > 4.0)]
                print(f"  [{tag}] SOLPS segments used: "
                      f"{len(solps_for_fixed)}/{len(solps_segments)}")
                map_solps_to_fw_bins(fb, solps_for_fixed)
                n_data = sum(1 for b in fb if b.data.get("Wtot", 0.0) != 0.0)
                print(f"  [{tag}] Bins with non-zero Wtot: {n_data}/{len(fb)}")

                out_coords = os.path.join(OUTPUT_DIR, f"mapped_{tag}_coords.dat")
                out_data   = os.path.join(OUTPUT_DIR, f"{tag}_data.dat")
                write_bin_coords(fb, out_coords)
                write_bin_data(fb, out_data)
                all_fixed_bins.extend(fb)
    else:
        print("\n── No fixed-bin files found in fixed_bins/. Skipping. ──")

    # ══════════════════════════════════════════════════════════════════════
    # COORDINATES TO BE BINNED  (Reactor_coordinates/coordinates_to_be_binned/)
    # Treatment: shadow-filter, subdivide, map, snap, merge.
    # ══════════════════════════════════════════════════════════════════════
    all_processed_bins: List[NewBin] = []
    all_div_groups: List[List[Tuple[float, float]]] = []

    if COORDS_TO_BIN_FILES:
        print(f"\n── Coordinates to be binned ({len(COORDS_TO_BIN_FILES)} file(s)) ──")
        for cb_path in COORDS_TO_BIN_FILES:
            tag = os.path.splitext(os.path.basename(cb_path))[0]
            print(f"\n  [{tag}] Parsing {os.path.basename(cb_path)} ...")
            div_groups = parse_divertor(cb_path)
            all_div_groups.extend(div_groups)
            print(f"  [{tag}] Groups: {len(div_groups)}  "
                  f"({sum(len(g) for g in div_groups)} points)")

            # 3. Filter SOLPS segments relevant to this region
            print(f"  [{tag}] Filtering SOLPS segments near region ...")
            solps_filtered = filter_solps_for_divertor(
                solps_segments, div_groups,
                max_dist=MAX_MAPPING_DIST, max_angle=MAX_ANGLE_DEG,
            )
            print(f"  [{tag}] Kept {len(solps_filtered)}/{len(solps_segments)} "
                  "SOLPS segments")

            # 4. Create initial bins (shadow-filtered)
            print(f"  [{tag}] Creating initial bins (shadow filter) ...")
            bins = create_initial_bins(div_groups, solps_filtered,
                                       max_shadow_dist=MAX_SHADOW_DIST)
            print(f"  [{tag}] Initial bins (illuminated edges): {len(bins)}")

            # 4.5 Snap discontinuity endpoints
            print(f"  [{tag}] Snapping discontinuity endpoints (1 cm) ...")
            bins, n_snaps, fixed_vertices = snap_discontinuity_endpoints_to_near_segments(
                bins,
                near_tol=0.01, connect_eps=1e-6,
                endpoint_margin=0.01, drop_if_shorter_than=1e-7,
                max_iters=6,
            )
            print(f"  [{tag}] Applied {n_snaps} snaps; bins now: {len(bins)}")
            if fixed_vertices:
                print(f"  [{tag}] Fixed vertices from snapping: {len(fixed_vertices)}")

            # 5. Subdivide long bins
            print(f"  [{tag}] Subdividing bins for mapping resolution ...")
            bins = subdivide_bins(bins, solps_filtered, fixed_vertices=fixed_vertices)
            print(f"  [{tag}] Bins after subdivision: {len(bins)}")

            # 6. Map SOLPS data onto bins
            print(f"  [{tag}] Mapping SOLPS data onto bins ...")
            map_solps_to_bins(bins, solps_segments,
                              max_dist=None, n_subsample=SUBSAMPLE_N)
            n_with_data = sum(1 for b in bins if b.data.get("Wtot", 0.0) != 0.0)
            print(f"  [{tag}] Bins with non-zero Wtot: {n_with_data}/{len(bins)}")

            # 6b. Save pre-merge bins
            pre_coords = os.path.join(OUTPUT_DIR, f"pre_merge_{tag}_coords.dat")
            pre_data   = os.path.join(OUTPUT_DIR, f"pre_merge_{tag}_data.dat")
            write_bin_coords(bins, pre_coords)
            write_bin_data(bins, pre_data)

            # 7. Merge similar adjacent bins
            print(f"  [{tag}] Merging similar adjacent bins ...")
            n_before = len(bins)
            bins = merge_similar_bins(
                bins,
                rtol=MERGE_RTOL,
                rtol_hot=MERGE_RTOL_HOT,
                hot_thresh=MERGE_HOT_THRESH,
                max_length=MAX_BIN_LENGTH,
                fixed_vertices=fixed_vertices,
            )
            print(f"  [{tag}] Bins after merging: {len(bins)} (was {n_before})")

            # Write region outputs
            out_coords = os.path.join(OUTPUT_DIR, f"{tag}_bins_coords.dat")
            out_data   = os.path.join(OUTPUT_DIR, f"{tag}_bins_data.dat")
            write_bin_coords(bins, out_coords)
            write_bin_data(bins, out_data)

            all_processed_bins.extend(bins)
    else:
        print("\n── No coordinate files found in coordinates_to_be_binned/. Skipping. ──")

    # ── Write combined outputs ──
    if all_processed_bins:
        write_bin_coords(all_processed_bins, OUT_COORDS)
        write_bin_data(all_processed_bins, OUT_DATA)

    # ── Conservation check (processed bins vs SOLPS in bounding box) ──
    pct_wtot = float('nan')
    pct_dflux = float('nan')
    if all_processed_bins:
        margin = 0.05
        bin_rs = [b.r_start for b in all_processed_bins] + [b.r_end for b in all_processed_bins]
        bin_zs = [b.z_start for b in all_processed_bins] + [b.z_end for b in all_processed_bins]
        r_lo, r_hi = min(bin_rs) - margin, max(bin_rs) + margin
        z_lo, z_hi = min(bin_zs) - margin, max(bin_zs) + margin

        solps_in_box = [s for s in solps_segments
                        if r_lo <= s.r_centre <= r_hi and z_lo <= s.z_centre <= z_hi]

        solps_wtot_integral = sum(
            seg.data.get("Wtot", 0.0) * seg_length(seg.r_start, seg.z_start, seg.r_end, seg.z_end)
            for seg in solps_in_box)
        solps_dflux_integral = sum(
            (seg.data.get("flxi_D", 0.0) + seg.data.get("flxa_D", 0.0))
            * seg_length(seg.r_start, seg.z_start, seg.r_end, seg.z_end)
            for seg in solps_in_box)
        bins_wtot_integral = sum(b.data.get("Wtot", 0.0) * b.arc_length for b in all_processed_bins)
        bins_dflux_integral = sum(
            (b.data.get("flxi_D", 0.0) + b.data.get("flxa_D", 0.0)) * b.arc_length
            for b in all_processed_bins)
        pct_wtot = 100.0 * bins_wtot_integral / solps_wtot_integral if solps_wtot_integral else float('nan')
        pct_dflux = 100.0 * bins_dflux_integral / solps_dflux_integral if solps_dflux_integral else float('nan')
        print(f"\n  Conservation ({len(solps_in_box)} SOLPS segs in bounding box):")
        print(f"    Wtot   – SOLPS: {solps_wtot_integral:.4e} W   Bins: {bins_wtot_integral:.4e} W   → {pct_wtot:.1f}%")
        print(f"    D flux – SOLPS: {solps_dflux_integral:.4e} 1/s  Bins: {bins_dflux_integral:.4e} 1/s  → {pct_dflux:.1f}%")

    # ── Summary ──
    print("\n" + "─" * 70)
    print(f"  Fixed bins   : {len(all_fixed_bins)}")
    print(f"  Processed bins: {len(all_processed_bins)}")
    if all_processed_bins:
        total_len = sum(b.arc_length for b in all_processed_bins)
        print(f"  Total arc-length (processed): {total_len:.4f} m")
        print(f"  Average bin length: {total_len / len(all_processed_bins):.4f} m")
        print(f"  Min bin length: {min(b.arc_length for b in all_processed_bins):.4f} m")
        print(f"  Max bin length: {max(b.arc_length for b in all_processed_bins):.4f} m")
    print("─" * 70)

    # ── Binned Flux CSV (combined fixed + processed) ──
    print(f"\n[CSV] Writing {os.path.basename(OUT_BINNED_FLUX_CSV)} ...")
    write_binned_flux_csv(all_processed_bins, all_fixed_bins, OUT_BINNED_FLUX_CSV)

    # ── Profile plots ──
    print("\nGenerating profile plots ...")
    if all_fixed_bins and all_processed_bins:
        plot_profiles(all_processed_bins, all_fixed_bins, OUTPUT_DIR,
                      solps_all=solps_segments)
    else:
        print("  Skipped profile plots (need both fixed and processed bins).")

    # ── Diagnostic plot ──
    print("\nGenerating diagnostic plot ...")
    solps_filtered_all = filter_solps_for_divertor(
        solps_segments, all_div_groups,
        max_dist=MAX_MAPPING_DIST, max_angle=MAX_ANGLE_DEG,
    ) if all_div_groups else []
    plot_all(solps_segments, solps_filtered_all, all_div_groups,
             all_processed_bins, OUT_PLOT,
             conservation=(pct_wtot, pct_dflux),
             fw_bins=all_fixed_bins if all_fixed_bins else None)

    print("\nDone.  All outputs in:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
