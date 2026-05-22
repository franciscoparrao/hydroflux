"""Extract the daily streamflow series of DGA station Rio Huasco En
Santa Juana (codigo 3820003) from the CR2 qflxDaily 2020 archive,
write it as a clean parquet file, and identify candidate flood
events for Track A calibration of the hydroflux 2D solver.

The CR2 archive layout (cr2_qflxDaily_2020):
- One row per metadata field (lines 1-15: codigo_estacion, ...,
  inicio_automatica).
- One row per date from 1900-01-01 onwards (lines 16+).
- One column per station (812 columns total).

Santa Juana sits at column 120 with code 03820003. Its observation
window is 1928-02-01 to 2019-07-31 (19 860 daily observations, 92
years — the longest gauge record in the Huasco basin).

Reproducir:
    python3 extract.py

Output:
    output/santa_juana_qflx.parquet  — clean daily series
    output/events_candidate.csv       — top events for calibration
    Console: summary statistics + flagged events.
"""

from __future__ import annotations

import zipfile
from pathlib import Path

import pandas as pd

CR2_ZIP = Path(
    "/home/franciscoparrao/proyectos/marea_roja/data/external/dga/cr2_qflxDaily_2019.zip"
)
STATION_CODE = "03820003"   # Rio Huasco En Santa Juana
STATION_NAME = "Rio Huasco En Santa Juana"
SUBDIR = "cr2_qflxDaily_2020"      # directory inside the zip
OUT_DIR = Path(__file__).parent / "output"


def extract_station_series(zip_path: Path, station_code: str) -> pd.DataFrame:
    """Read the date column + the single station column from the CR2
    archive (without extracting the full 218 MB)."""
    with zipfile.ZipFile(zip_path) as zf:
        data_name = f"{SUBDIR}/{SUBDIR}.txt"
        with zf.open(data_name) as f:
            df = pd.read_csv(
                f,
                usecols=["codigo_estacion", station_code],
                dtype={"codigo_estacion": "string", station_code: "string"},
                low_memory=False,
            )
    # First 15 rows are metadata (codigo_estacion column holds field
    # names like "institucion", "altura", ...). The rest hold dates.
    metadata = df.iloc[:15].set_index("codigo_estacion").to_dict()[station_code]
    body = df.iloc[15:].copy()
    body = body.rename(columns={"codigo_estacion": "date", station_code: "qflx_m3s"})
    body["date"] = pd.to_datetime(body["date"], format="%Y-%m-%d", errors="coerce")
    # CR2 uses both empty cells AND the sentinel -9999 for missing
    # observations (legacy DGA format). Both must collapse to NaN
    # for the downstream summary statistics to be meaningful.
    body["qflx_m3s"] = pd.to_numeric(body["qflx_m3s"], errors="coerce")
    body.loc[body["qflx_m3s"] <= -9000, "qflx_m3s"] = pd.NA
    body = body.dropna(subset=["date"]).reset_index(drop=True)
    body.attrs["metadata"] = metadata
    body.attrs["station_code"] = station_code
    body.attrs["station_name"] = STATION_NAME
    return body


def summarise(df: pd.DataFrame) -> None:
    meta = df.attrs["metadata"]
    valid = df.dropna(subset=["qflx_m3s"])
    print(f"\n=== {STATION_NAME} (codigo {STATION_CODE}) ===")
    print(f"Coordinates: lat={meta.get('latitud')}, lon={meta.get('longitud')}")
    print(f"Altitude:    {meta.get('altura')} m")
    print(f"Sub-basin:   {meta.get('nombre_sub_cuenca')} (code {meta.get('codigo_sub_cuenca')})")
    print(f"Source:      {meta.get('institucion')} via {meta.get('fuente')}")
    print()
    print(f"Date range:      {df['date'].min().date()} → {df['date'].max().date()}")
    print(f"Total rows:      {len(df):,}")
    print(f"Valid (non-NaN): {len(valid):,}  ({100 * len(valid) / len(df):.1f}%)")
    print()
    if len(valid):
        q = valid["qflx_m3s"]
        print(f"q [m³/s]   min  : {q.min():.3f}")
        print(f"           p05  : {q.quantile(0.05):.3f}")
        print(f"           p50  : {q.quantile(0.50):.3f}")
        print(f"           mean : {q.mean():.3f}")
        print(f"           p95  : {q.quantile(0.95):.3f}")
        print(f"           p99  : {q.quantile(0.99):.3f}")
        print(f"           max  : {q.max():.3f}")


def identify_events(
    df: pd.DataFrame,
    n_top: int = 10,
    min_separation_days: int = 60,
) -> pd.DataFrame:
    """Pick the `n_top` highest-flow daily peaks with at least
    `min_separation_days` between each, so we get distinct events
    rather than several days of one event."""
    valid = df.dropna(subset=["qflx_m3s"]).reset_index(drop=True)
    valid = valid.sort_values("qflx_m3s", ascending=False).reset_index(drop=True)

    selected: list[pd.Series] = []
    for _, row in valid.iterrows():
        if all(
            abs((row["date"] - s["date"]).days) >= min_separation_days
            for s in selected
        ):
            selected.append(row)
            if len(selected) >= n_top:
                break

    out = pd.DataFrame(selected).reset_index(drop=True)
    out["year"] = out["date"].dt.year
    out["month"] = out["date"].dt.month
    return out[["date", "year", "month", "qflx_m3s"]]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not CR2_ZIP.exists():
        raise SystemExit(f"CR2 archive not found at {CR2_ZIP}")

    print(f"Reading {CR2_ZIP} ...")
    df = extract_station_series(CR2_ZIP, STATION_CODE)

    parquet_path = OUT_DIR / "santa_juana_qflx.parquet"
    df.to_parquet(parquet_path, index=False)
    print(f"Wrote daily series → {parquet_path}")

    summarise(df)

    events = identify_events(df)
    events_path = OUT_DIR / "events_candidate.csv"
    events.to_csv(events_path, index=False)

    print(f"\n=== Top {len(events)} peak-flow candidate events (min 60 d apart) ===")
    print(events.to_string(index=False))
    print(f"\nWrote candidate events → {events_path}")


if __name__ == "__main__":
    main()
