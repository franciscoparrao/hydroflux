"""Extract a DEM-derived longitudinal profile of the Río Huasco at
Santa Juana for use as input to the hydroflux 1D solver.

Procedure
---------
1. Reproject the gauge location (lon, lat WGS84) into the DEM CRS
   (UTM 19S, EPSG:32719).
2. Snap the gauge to the highest-flow-accumulation cell within a
   small neighbourhood — this puts us ON the main stem, not on a
   tributary or on the bank.
3. Walk downstream along the D8 flow direction for N cells (≈ N·30 m
   along channel). At each step, record the cell's bed elevation
   (from the pit-filled DEM, which avoids spurious local minima
   that would break the monotone descending assumption of the 1D
   solver).
4. Convert the (x, y, z) trajectory into a (distance, elevation)
   profile using cumulative arc length between consecutive cell
   centres.
5. Compute the per-cell slope as a finite difference on the
   profile. Also fit a single best-fit linear slope across the
   whole reach (useful as a sanity check vs the synthetic 0.005
   used in the previous demo).

Output
------
- `huasco_longitudinal_profile.csv`: columns
   `distance_m, elevation_m, slope` (slope is forward-difference; NaN at the
   last sample).
- `huasco_profile_summary.json`: aggregate stats (total length, mean
   slope, fitted slope, elevation drop) — handy for the Rust demo's
   header.

Reproducir:
    python3 extract_longitudinal_profile.py [--n-cells 50]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import transform as warp_transform

GAUGE_LON = -70.6464
GAUGE_LAT = -28.6719
GAUGE_CRS = "EPSG:4326"

FACTORS_DIR = Path(
    "/home/franciscoparrao/proyectos/postdoc/papers/paper1_susceptibilidad/factors/06_rio_huasco"
)
DEM_PATH = FACTORS_DIR / "hydrology" / "filled.tif"
FLOW_DIR_PATH = FACTORS_DIR / "hydrology" / "flow_direction_d8.tif"
FLOW_ACC_PATH = FACTORS_DIR / "hydrology" / "flow_accumulation.tif"
HAND_PATH = FACTORS_DIR / "hydrology" / "hand.tif"

# HAND threshold (metres) used to define "in-channel" for the
# connected-perpendicular width estimate. 0.5 m corresponds roughly
# to the active channel at low-to-moderate flow; higher thresholds
# include the floodplain. Picked conservatively because the
# pit-filled DEM creates artificially flat HAND=0 pools that would
# inflate the perpendicular extent if we walked through them
# disconnectedly.
HAND_CHANNEL_THRESHOLD_M = 0.5
MAX_WIDTH_HALF_CELLS = 30

# D8 flow-direction encoding used by the susceptibility pipeline's
# WhiteboxTools-via-TauDEM output: integers 1-8 where 1 = E, then
# counter-clockwise around the compass:
#   4  3  2
#   5     1
#   6  7  8
# 0 = no defined flow (sink / nodata). Empirically verified for
# (906, 1791) where elevation 490.50 m has W=488.85 (lowest neighbour)
# and the raster value at that cell is 5.
D8_OFFSETS = {
    1: (0, 1),    # E
    2: (-1, 1),   # NE
    3: (-1, 0),   # N
    4: (-1, -1),  # NW
    5: (0, -1),   # W
    6: (1, -1),   # SW
    7: (1, 0),    # S
    8: (1, 1),    # SE
}


def find_gauge_pixel(dem):
    """Project the gauge to the DEM CRS and return the (row, col) of
    the containing cell."""
    xs, ys = warp_transform(GAUGE_CRS, str(dem.crs), [GAUGE_LON], [GAUGE_LAT])
    row, col = dem.index(xs[0], ys[0])
    return int(row), int(col)


def snap_to_main_stem(flow_acc, row, col, window=5):
    """Move the (row, col) to the highest-accumulation cell within a
    `(2·window+1)²` window. Keeps us on the actual main stem if the
    gauge coordinate lands on the bank."""
    h, w = flow_acc.shape
    r_lo, r_hi = max(0, row - window), min(h, row + window + 1)
    c_lo, c_hi = max(0, col - window), min(w, col + window + 1)
    patch = flow_acc[r_lo:r_hi, c_lo:c_hi]
    flat = int(patch.argmax())
    dr, dc = np.unravel_index(flat, patch.shape)
    return r_lo + int(dr), c_lo + int(dc)


def channel_width_at(
    hand,
    flow_dir,
    row,
    col,
    pixel_size_m=30.0,
    threshold=HAND_CHANNEL_THRESHOLD_M,
    max_half_cells=MAX_WIDTH_HALF_CELLS,
):
    """Connected-perpendicular channel width at (row, col).

    Compute the perpendicular direction to the cell's D8 flow vector,
    then walk OUT in both perpendicular senses until the first cell
    with HAND > threshold (or the raster edge / nodata). The width is
    the total number of consecutive in-channel cells × pixel step.

    Connected walk: stops at the first out-of-channel cell, so it
    cannot bleed into disconnected flat HAND=0 pools elsewhere in
    the basin (a problem the pit-filled DEM creates).

    For diagonal perpendiculars (flow direction NE/SE/SW/NW), the
    perpendicular step length is `pixel · √2`; we account for that
    in the width sum.
    """
    d = int(flow_dir[row, col])
    if d not in D8_OFFSETS:
        return float("nan")
    dr_f, dc_f = D8_OFFSETS[d]
    # Two perpendiculars: (-dc, dr) and (dc, -dr).
    pdr, pdc = -dc_f, dr_f
    step_len = pixel_size_m * (1.0 if (abs(pdr) + abs(pdc) == 1) else 2**0.5)
    ncells = 1  # center cell

    def walk(sign):
        r, c = row, col
        out = 0
        for _ in range(max_half_cells):
            r, c = r + sign * pdr, c + sign * pdc
            if not (0 <= r < hand.shape[0] and 0 <= c < hand.shape[1]):
                break
            h = hand[r, c]
            if np.isnan(h) or h > threshold:
                break
            out += 1
        return out

    ncells += walk(+1) + walk(-1)
    return ncells * step_len


def walk_downstream(flow_dir, dem, start_row, start_col, n_cells):
    """Walk along D8 flow direction for at most `n_cells` steps,
    recording (row, col, x, y, z) at each step. Stops early if the
    flow direction is invalid (edge of basin) or if the path turns
    upstream (z increases — defensive against degenerate cells)."""
    row, col = start_row, start_col
    trajectory = []
    h, w = flow_dir.shape
    transform = dem.transform
    z_prev = dem.read(1)[row, col]
    for _ in range(n_cells):
        x, y = transform * (col + 0.5, row + 0.5)
        z = dem.read(1)[row, col]
        trajectory.append((row, col, float(x), float(y), float(z)))
        d = int(flow_dir[row, col])
        if d not in D8_OFFSETS:
            break
        dr, dc = D8_OFFSETS[d]
        nrow, ncol = row + dr, col + dc
        if not (0 <= nrow < h and 0 <= ncol < w):
            break
        z_next = dem.read(1)[nrow, ncol]
        if z_next > z_prev + 0.5:
            # Defensive: D8 sometimes loops on flat patches. The
            # pit-filled DEM should avoid this, but bail if it
            # happens anyway.
            break
        row, col, z_prev = nrow, ncol, z_next
    return trajectory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--n-cells", type=int, default=50, help="Number of D8 steps downstream (default 50 ≈ 1500 m)"
    )
    parser.add_argument(
        "--snap-window",
        type=int,
        default=150,
        help=(
            "Half-window (cells) for main-stem snap (default 150 ≈ 4.5 km). "
            "Larger windows tolerate gauge-coordinate / DEM datum offsets "
            "(common in Chile where DGA coordinates may be PSAD56 vs the "
            "WGS84 DEM); chosen large enough to find the Huasco main stem "
            "for the Santa Juana gauge."
        ),
    )
    parser.add_argument(
        "--min-acc",
        type=float,
        default=5.0e6,
        help=(
            "Minimum flow-accumulation cells the snapped pixel must have. "
            "Santa Juana sits below ~7700 km² of contributing area on the "
            "main Huasco, so the snapped pixel should report ~8M cells. "
            "Bail if we end up below this — usually a sign the gauge "
            "coordinate is wrong or the DEM extent does not actually "
            "contain Santa Juana."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).parent / "output",
        help="Where to write the CSV + JSON",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    with rasterio.open(DEM_PATH) as dem, rasterio.open(FLOW_DIR_PATH) as fdir, rasterio.open(
        FLOW_ACC_PATH
    ) as facc, rasterio.open(HAND_PATH) as hand_src:
        # Find and snap the gauge.
        row0, col0 = find_gauge_pixel(dem)
        print(f"Gauge raw pixel: ({row0}, {col0})")
        flow_acc = facc.read(1)
        row, col = snap_to_main_stem(flow_acc, row0, col0, window=args.snap_window)
        snapped_acc = float(flow_acc[row, col])
        print(
            f"Snapped to main stem: ({row}, {col}); "
            f"flow_acc raw = {flow_acc[row0, col0]:.0f}, snapped = {snapped_acc:.0f}"
        )
        if snapped_acc < args.min_acc:
            raise SystemExit(
                f"Snapped pixel has acc={snapped_acc:.0f} < --min-acc={args.min_acc:.0f}. "
                "Probable cause: gauge coordinate is off (try widening --snap-window or "
                "double-checking the DGA coordinate datum)."
            )

        # Walk.
        flow_dir = fdir.read(1)
        hand = hand_src.read(1)
        traj = walk_downstream(flow_dir, dem, row, col, args.n_cells)

        # Per-cell channel width via connected-perpendicular HAND walk.
        widths = np.array(
            [channel_width_at(hand, flow_dir, r, c) for (r, c, _, _, _) in traj]
        )

    if len(traj) < 2:
        raise SystemExit(
            f"Trajectory has only {len(traj)} sample(s) — D8 walk stopped immediately. "
            "Likely the snapped pixel has invalid flow direction or is at the basin outlet."
        )

    # Build longitudinal profile.
    xs = np.array([p[2] for p in traj])
    ys = np.array([p[3] for p in traj])
    zs = np.array([p[4] for p in traj])
    dx_arc = np.hypot(np.diff(xs), np.diff(ys))
    distance = np.concatenate(([0.0], np.cumsum(dx_arc)))

    # Per-segment slope (forward difference). Last cell gets NaN.
    slope = np.full_like(zs, np.nan)
    slope[:-1] = -np.diff(zs) / dx_arc  # positive when bed descends in +x

    # CSV — now includes width.
    csv_path = args.output_dir / "huasco_longitudinal_profile.csv"
    with csv_path.open("w") as f:
        f.write("distance_m,elevation_m,slope,width_m\n")
        for d, z, s, w in zip(distance, zs, slope, widths):
            sf = "" if np.isnan(s) else f"{s:.6f}"
            wf = "" if np.isnan(w) else f"{w:.2f}"
            f.write(f"{d:.3f},{z:.4f},{sf},{wf}\n")
    print(f"Wrote {csv_path} ({len(distance)} samples, with width)")

    # Aggregate summary.
    total_length = float(distance[-1])
    elev_drop = float(zs[0] - zs[-1])
    mean_slope = elev_drop / total_length if total_length > 0 else 0.0
    valid_slopes = slope[:-1][~np.isnan(slope[:-1])]
    seg_slope_median = float(np.median(valid_slopes)) if len(valid_slopes) > 0 else 0.0
    # Best-fit linear slope through the (distance, elevation) cloud.
    poly = np.polyfit(distance, zs, 1)
    fitted_slope = float(-poly[0])  # negative gradient → positive descent rate

    valid_widths = widths[~np.isnan(widths)]
    width_median = float(np.median(valid_widths)) if len(valid_widths) > 0 else float("nan")
    width_mean = float(np.mean(valid_widths)) if len(valid_widths) > 0 else float("nan")
    width_p25 = float(np.percentile(valid_widths, 25)) if len(valid_widths) > 0 else float("nan")
    width_p75 = float(np.percentile(valid_widths, 75)) if len(valid_widths) > 0 else float("nan")

    summary = {
        "gauge_lon": GAUGE_LON,
        "gauge_lat": GAUGE_LAT,
        "n_cells_walked": len(traj),
        "total_length_m": total_length,
        "elevation_drop_m": elev_drop,
        "mean_slope": mean_slope,
        "segment_slope_median": seg_slope_median,
        "fitted_slope_linear": fitted_slope,
        "elevation_start_m": float(zs[0]),
        "elevation_end_m": float(zs[-1]),
        "hand_threshold_m": HAND_CHANNEL_THRESHOLD_M,
        "channel_width_median_m": width_median,
        "channel_width_mean_m": width_mean,
        "channel_width_p25_m": width_p25,
        "channel_width_p75_m": width_p75,
        "note_resolution_limit": (
            "DEM 30 m: the narrowest resolvable channel is 1 pixel = 30 m. "
            "If the true active channel at the gauge is narrower (e.g., 5-15 m), "
            "the DEM-derived width will overestimate it. A higher-resolution DEM "
            "(LiDAR / Pleiades 0.5 m) is the only fix at this stage."
        ),
    }
    json_path = args.output_dir / "huasco_profile_summary.json"
    json_path.write_text(json.dumps(summary, indent=2))
    print(f"Wrote {json_path}")
    print(json.dumps(summary, indent=2))

    # Regridded uniform-dx profile for direct consumption by the 1D
    # solver. The hydroflux Lax-Friedrichs scheme assumes uniform dx;
    # the raw D8 walk alternates 30 m (cardinal) and 42 m (diagonal)
    # steps. Linear interpolation onto a uniform grid handles both
    # the dx irregularity and any small per-cell elevation noise
    # (DEM filled artefacts) without changing the integrated bed
    # drop or the gross slope.
    n_uniform = 60
    distance_uniform = np.linspace(0.0, total_length, n_uniform)
    elev_uniform = np.interp(distance_uniform, distance, zs)
    rust_const = ", ".join(f"{z:.4f}" for z in elev_uniform)
    rust_snippet = (
        f"// {n_uniform} samples uniformly spaced over {total_length:.1f} m "
        f"(dx = {total_length / (n_uniform - 1):.4f} m).\n"
        f"// Mean slope {mean_slope:.6f}, fitted linear {fitted_slope:.6f}, "
        f"drop {elev_drop:.3f} m.\n"
        f"// Generated by examples/huasco_channel/extract_longitudinal_profile.py.\n"
        f"const HUASCO_BED_M: [f64; {n_uniform}] = [\n    {rust_const},\n];\n"
    )
    rust_path = args.output_dir / "huasco_bed.rs.snippet"
    rust_path.write_text(rust_snippet)
    print(f"Wrote {rust_path} (paste into Rust example)")


if __name__ == "__main__":
    main()
