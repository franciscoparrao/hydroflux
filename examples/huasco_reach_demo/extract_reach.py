"""Extract a 1D bed profile along a tributary of Río Huasco.

Contrasts with the Maule demo: this basin is semiarid Andean (BNA #06),
with steeper terrain (peaks > 6000 m), episodic flows, and a much drier
hydrology. The same D8 tracing logic is used; only the filter ranges
shift to target the proper terrain regime.

Same pipeline shape as `maule_reach_demo/extract_reach.py`.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import Affine

DATA = (
    Path.home()
    / "proyectos/postdoc/papers/paper1_susceptibilidad/factors/06_rio_huasco"
)
OUT_DIR = Path(__file__).parent / "output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TARGET_LEN_M = 10_000.0
# Huasco-specific filter: higher elevation band (Andean proper, headwaters
# of a semiarid basin where main channels sit above 1500 m), smaller
# accumulation bracket since episodic flows produce smaller catchments.
ELEV_RANGE = (1500.0, 3500.0)
ACCUM_RANGE = (20_000, 200_000)

# TauDEM D8 (verified shared with the Maule data).
D8_DR = {1: 0, 2: -1, 3: -1, 4: -1, 5: 0, 6: 1, 7: 1, 8: 1}
D8_DC = {1: 1, 2: 1, 3: 0, 4: -1, 5: -1, 6: -1, 7: 0, 8: 1}


def main() -> None:
    with rasterio.open(DATA / "hydrology" / "filled.tif") as src:
        dem = src.read(1).astype(np.float64)
        nd = src.nodata
        if nd is not None:
            dem[dem == nd] = np.nan
        crs = src.crs
        dem_transform = src.transform

    with rasterio.open(DATA / "hydrology" / "flow_direction_d8.tif") as src:
        fd = src.read(1).astype(np.int32)

    with rasterio.open(DATA / "hydrology" / "flow_accumulation.tif") as src:
        fa = src.read(1).astype(np.float64)

    rows, cols = fd.shape
    pixel_size = abs(dem_transform.a)
    print(f"DEM: {rows}×{cols} cells, pixel size {pixel_size} m")

    candidates = (
        (fd >= 1)
        & (fd <= 8)
        & ~np.isnan(dem)
        & (dem >= ELEV_RANGE[0])
        & (dem <= ELEV_RANGE[1])
        & (fa >= ACCUM_RANGE[0])
        & (fa <= ACCUM_RANGE[1])
    )
    if not candidates.any():
        raise SystemExit("no cells matched the candidate filter — relax ranges")
    masked = np.where(candidates, fa, -1.0)
    r0, c0 = np.unravel_index(int(np.argmax(masked)), masked.shape)
    print(
        f"Start: row={r0}, col={c0}, accum={fa[r0, c0]:.0f} cells "
        f"(~{fa[r0, c0] * pixel_size ** 2 / 1e6:.1f} km²)"
    )
    print(f"Start elevation: {dem[r0, c0]:.1f} m")

    path_rc: list[tuple[int, int]] = [(r0, c0)]
    visited = {(r0, c0)}
    total_len = 0.0
    diag_len = pixel_size * np.sqrt(2.0)
    prev_acc = fa[r0, c0]
    r, c = r0, c0
    while total_len < TARGET_LEN_M:
        d = int(fd[r, c])
        if d not in D8_DR:
            print(f"Stopped: invalid D8 code {d} at ({r}, {c})")
            break
        dr, dc = D8_DR[d], D8_DC[d]
        nr, nc = r + dr, c + dc
        if not (0 <= nr < rows and 0 <= nc < cols):
            print(f"Stopped: out of bounds at ({nr}, {nc})")
            break
        if (nr, nc) in visited:
            print(f"Stopped: cycle at ({nr}, {nc})")
            break
        if np.isnan(dem[nr, nc]):
            print(f"Stopped: NaN bed at ({nr}, {nc})")
            break
        step = diag_len if abs(dr) + abs(dc) == 2 else pixel_size
        total_len += step
        path_rc.append((nr, nc))
        visited.add((nr, nc))
        if fa[nr, nc] < prev_acc * 0.5:
            print(
                f"Warning: accumulation dropped sharply at ({nr}, {nc}): "
                f"{prev_acc:.0f} → {fa[nr, nc]:.0f} (possible bad D8 reading)"
            )
        prev_acc = fa[nr, nc]
        r, c = nr, nc

    n_cells = len(path_rc)
    print(f"Reach: {n_cells} cells, total length {total_len:.0f} m")

    bed = np.array([dem[r, c] for r, c in path_rc], dtype=np.float32)
    print(
        f"Bed: top {bed[0]:.1f} m, bottom {bed[-1]:.1f} m, "
        f"drop {bed[0] - bed[-1]:.1f} m, mean slope "
        f"{(bed[0] - bed[-1]) / total_len:.4f}"
    )

    dx = total_len / (n_cells - 1)
    print(f"Effective dx = {dx:.2f} m")

    csv_path = OUT_DIR / "centerline.csv"
    xs = np.cumsum([0.0] + [
        diag_len if abs(path_rc[i + 1][0] - path_rc[i][0])
        + abs(path_rc[i + 1][1] - path_rc[i][1]) == 2
        else pixel_size
        for i in range(n_cells - 1)
    ])
    with open(csv_path, "w") as f:
        f.write("station_m,row,col,bed_m\n")
        for x, (r, c), b in zip(xs, path_rc, bed):
            f.write(f"{x:.3f},{r},{c},{float(b):.4f}\n")
    print(f"Wrote {csv_path}")

    out_path = OUT_DIR / "bed.tif"
    transform = Affine(dx, 0.0, 0.0, 0.0, -dx, 0.0)
    with rasterio.open(
        out_path,
        "w",
        driver="GTiff",
        height=1,
        width=n_cells,
        count=1,
        dtype=rasterio.float32,
        crs=crs,
        transform=transform,
    ) as dst:
        dst.write(bed.reshape(1, -1), 1)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
