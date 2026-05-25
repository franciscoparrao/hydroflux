"""Extract a 6 km × 2 km subset of the Huasco DEM centred on the
Santa Juana gauge for the solver-2d Phase 2 first real-data
simulation. Also emits a co-registered flow-accumulation raster so
the Rust example can identify channel cells at the East inflow edge.

# Window definition

Snap point (from the longitudinal-profile script): (row, col) =
(906, 1791) in the full DEM, corresponding to the Santa Juana gauge
re-projected to UTM 19S and snapped to the highest-flow-accumulation
cell within a 4.5 km window.

Window: 200 cols × 67 rows = 6000 m × 2010 m around the snap point,
keeping the main Huasco stem running roughly E→W through the
centre. Cols [1691, 1891), rows [873, 940), giving exactly
200 × 67 cells.

Output: `huasco_subset_dem.tif` and `huasco_subset_acc.tif` in the
same window with the same GeoTransform (so the Rust loader can use
them interchangeably).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window

FACTORS_DIR = Path(
    "/home/franciscoparrao/proyectos/postdoc/papers/paper1_susceptibilidad/factors/06_rio_huasco"
)
DEM_PATH = FACTORS_DIR / "hydrology" / "filled.tif"
FLOW_ACC_PATH = FACTORS_DIR / "hydrology" / "flow_accumulation.tif"

# Snap point identified by examples/huasco_channel/extract_longitudinal_profile.py.
# Hard-coded here so this script is self-contained.
SNAP_ROW = 906
SNAP_COL = 1791

# Window dimensions (cells).
#
# The Huasco main stem at the Santa Juana gauge runs **roughly
# south-to-north** within a 6 km window (verified empirically:
# accumulation grows by ~10⁴ cells from S to N inside the subset,
# meaning downstream is N). To get clean inflow/outflow on the short
# edges of the window, the rectangle is oriented **portrait**: tall
# along the river direction (N-S) and narrow across (E-W).
#
# 67 cols × 200 rows = 2010 m (E-W) × 6000 m (N-S) at 30 m DEM.
WIN_COLS = 67
WIN_ROWS = 200


def extract_window(src_path: Path, dst_path: Path, window: Window) -> None:
    """Read `src_path`, slice `window`, write to `dst_path` updating
    the GeoTransform to the new origin."""
    with rasterio.open(src_path) as src:
        data = src.read(1, window=window)
        new_transform = src.window_transform(window)
        profile = src.profile.copy()
        profile.update(
            height=int(window.height),
            width=int(window.width),
            transform=new_transform,
            count=1,
            # Cast to float32 — what the SurtGIS native GeoTIFF writer
            # produces and what solver-2d::io reads cleanly.
            dtype="float32",
            compress="lzw",
        )
        with rasterio.open(dst_path, "w", **profile) as dst:
            dst.write(data.astype("float32"), 1)
    print(f"Wrote {dst_path} (shape {window.height}×{window.width}, "
          f"origin = {new_transform.c:.1f}, {new_transform.f:.1f})")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--center-row", type=int, default=SNAP_ROW,
        help=f"Window center row in full DEM (default {SNAP_ROW}).",
    )
    parser.add_argument(
        "--center-col", type=int, default=SNAP_COL,
        help=f"Window center col in full DEM (default {SNAP_COL}).",
    )
    parser.add_argument(
        "--rows", type=int, default=WIN_ROWS,
        help=f"Window height in cells (default {WIN_ROWS}).",
    )
    parser.add_argument(
        "--cols", type=int, default=WIN_COLS,
        help=f"Window width in cells (default {WIN_COLS}).",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path(__file__).parent / "output",
        help="Where to write subset rasters.",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    row_off = args.center_row - args.rows // 2
    col_off = args.center_col - args.cols // 2
    window = Window(col_off=col_off, row_off=row_off,
                    width=args.cols, height=args.rows)
    print(f"Subset window: rows [{row_off}, {row_off + args.rows}), "
          f"cols [{col_off}, {col_off + args.cols}) "
          f"= {args.cols * 30}m × {args.rows * 30}m")

    extract_window(DEM_PATH,      args.output_dir / "huasco_subset_dem.tif", window)
    extract_window(FLOW_ACC_PATH, args.output_dir / "huasco_subset_acc.tif", window)

    # Identify the inflow + outflow cells by scanning all 4 edges
    # for channel cells (acc > threshold). Outflow = boundary cell
    # with highest accumulation (most downstream). Inflow = boundary
    # cells on the OPPOSITE side from outflow (the river enters
    # somewhere upstream of the outflow).
    with rasterio.open(args.output_dir / "huasco_subset_acc.tif") as src:
        acc = src.read(1)
    threshold = 1_000_000.0  # ~900 km² catchment → main Huasco stem
    channel = acc > threshold
    print(
        f"\nChannel cells (acc > {threshold:.0e}): {channel.sum()} total in subset.\n"
    )

    n_rows, n_cols = acc.shape
    edge_cells = []  # list of (edge_name, row, col, acc)
    for c in range(n_cols):
        if channel[0, c]:  # N edge (row 0 = top of image)
            edge_cells.append(("N", 0, c, acc[0, c]))
        if channel[n_rows - 1, c]:  # S edge (last row = bottom)
            edge_cells.append(("S", n_rows - 1, c, acc[n_rows - 1, c]))
    for r in range(n_rows):
        if channel[r, 0]:  # W edge (col 0)
            edge_cells.append(("W", r, 0, acc[r, 0]))
        if channel[r, n_cols - 1]:  # E edge (last col)
            edge_cells.append(("E", r, n_cols - 1, acc[r, n_cols - 1]))

    if not edge_cells:
        print(
            "ERROR: no channel cells on ANY boundary edge. The window is "
            "fully internal. Shift the window or use a smaller threshold."
        )
        return

    edge_cells.sort(key=lambda x: x[3], reverse=True)
    outflow = edge_cells[0]
    inflow_candidates = [c for c in edge_cells if c[0] != outflow[0]]

    print(f"OUTFLOW (highest-acc boundary cell): {outflow}")
    print(f"  → suggest setting BC '{outflow[0]}' edge as Transmissive\n")
    print(f"INFLOW candidates (other-edge channel cells): {len(inflow_candidates)}")
    for c in inflow_candidates:
        print(f"  ({c[1]:>3}, {c[2]:>3}) on edge {c[0]}: acc = {c[3]:.0f}")
    inflow_edges = sorted(set(c[0] for c in inflow_candidates))
    print(
        f"  → suggest injecting Q via PointSource at these cells "
        f"(inflow edges: {inflow_edges})"
    )
    print(f"  → suggest setting all other edges as Wall")


if __name__ == "__main__":
    main()
