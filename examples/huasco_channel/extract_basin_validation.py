"""Extract Q observations from all Huasco-basin DGA stations for the
2017 Atacama event window, for validation-split analysis of Track A
calibration at Santa Juana.

# Background

Track A application calibrates Manning n on a 1.8-km Huasco reach
just below the Santa Juana gauge (DGA code 3820003, lat -28.6719,
elevation 575 m). For physical defensibility, the calibration
should ideally be validated against an INDEPENDENT gauge — either
upstream (input forcing) or downstream (routing target).

This script extracts the 2017-02-20 → 2017-03-12 window for all
nearby DGA stations in the Huasco basin to assess what validation
is actually possible.

# Empirical finding (2026-05-24)

For the Atacama 2017 event window at Santa Juana (peak 38.9 m³/s,
mean 31.85 m³/s), the upstream tributary stations show only
BASEFLOW:

  Station                            Elev    Mean Q   Peak Q   Event?
  --------------------------------------------------------------------
  Santa Juana                        575 m   31.85    38.90    YES
  Carmen Ramadillas                  825 m    3.71     4.62    no
  Carmen Pte La Majada              1075 m    4.76     5.37    no
  Carmen El Corral                  2000 m    4.64     4.83    no
  Tránsito Angostura Pinte          1000 m    1.51     2.16    no
  Conay Las Lozas (Tránsito basin)            4.73     5.26    no
  Huasco Pte Nicolasa (downstream)   167 m   10.43    11.90    attenuated
  Tránsito Antes Junta Carmen        812 m    -        -       no 2017 data

Total upstream baseflow ≈ 18-19 m³/s, but Santa Juana measures
31-38 m³/s. The missing ~15-20 m³/s during the event came from
sub-basins between the Tránsito/Carmen confluence and Santa Juana
that have NO DGA gauges. Atacama 2017 was a localised event, not a
propagated wave from headwaters.

Downstream attenuation is also striking: Pte Nicolasa (167 m,
~50 km downstream of Santa Juana) sees mean 10.43 m³/s versus
31.85 at Santa Juana — typical heavy irrigation abstraction in
the Huasco agricultural valley.

# Implications for Track A

A clean upstream/downstream validation split is NOT POSSIBLE for
the 2017 event with these gauges. The compound-section calibration
on Santa Juana alone (iter 6, `calibrate_manning_huasco_2017_compound`)
is the most we can defend with these data.

For full validation a future iteration would need either:
1. Sub-daily data (where intra-day variability gives spatial
   information beyond what daily means provide), OR
2. Sub-basin gauges between confluence and Santa Juana (not
   available), OR
3. A precipitation-runoff model to derive synthetic upstream
   forcing for the un-gauged sub-basins (paper 2+ scope).

Alternatively, pick a DIFFERENT historic event where Tránsito
Antes Junta (record 1927-2015) and Carmen had concurrent peaks.
Top candidates from the original events_candidate.csv:
  1984-07-11 (Q=107 m³/s, max histórico) — pre-2015, all stations
  1998-01-07 (Q=93.6, La Niña fuerte) — pre-2015
  1984-12-05 (Q=84.9) — pre-2015

Reproducir:
    python3 extract_basin_validation.py
"""

import io
import json
import zipfile
from pathlib import Path

import numpy as np

CR2_ZIP = Path("/home/franciscoparrao/proyectos/marea_roja/data/external/dga/cr2_qflxDaily_2019.zip")
DATA_ENTRY = "cr2_qflxDaily_2020/cr2_qflxDaily_2020.txt"

# All Huasco-basin DGA stations with reasonable record lengths.
# Codes carry a leading zero in the CR2 column header (e.g. "03820003").
STATIONS = {
    "03820003": ("Santa Juana", 575),
    "03825001": ("Huasco Pte Nicolasa", 167),
    "03815001": ("Carmen Ramadillas", 825),
    "03815002": ("Carmen Pte La Majada", 1075),
    "03814003": ("Carmen El Corral", 2000),
    "03806001": ("Tránsito Antes Junta", 812),
    "03804002": ("Tránsito Angostura Pinte", 1000),
    "03802001": ("Conay Las Lozas", 0),
}

WIN_START = "2017-02-20"
WIN_END = "2017-03-12"


def extract_window(zip_path: Path, entry: str, codes: dict, win_start: str, win_end: str):
    with zipfile.ZipFile(zip_path) as z, z.open(entry) as f:
        f_text = io.TextIOWrapper(f, encoding="latin-1")
        header = f_text.readline().strip().split(",")
        col_codes = header[1:]
        targets = {c: col_codes.index(c) + 1 for c in codes if c in col_codes}
        missing = set(codes) - set(targets)
        if missing:
            print(f"NOT in CR2 data: {[(c, codes[c][0]) for c in missing]}")
        for _ in range(14):
            f_text.readline()
        rows = []
        for line in f_text:
            cells = line.strip().split(",")
            d = cells[0]
            if d < win_start:
                continue
            if d > win_end:
                break
            row = {"date": d}
            for c, idx in targets.items():
                s = cells[idx].strip() if idx < len(cells) else ""
                row[c] = float("nan") if s in ("", "-9999", "-9999.0") else float(s)
            rows.append(row)
    return rows, targets


def main():
    rows, targets = extract_window(CR2_ZIP, DATA_ENTRY, STATIONS, WIN_START, WIN_END)
    print(f"\nAtacama 2017 event window: {WIN_START} → {WIN_END}  ({len(rows)} days)\n")
    print(f"{'Date':<12}" + "".join(f"{STATIONS[c][0][:18]:>20}" for c in targets))
    for row in rows:
        parts = [f"{row['date']:<12}"]
        for c in targets:
            v = row[c]
            parts.append(f"{v:>20.2f}" if not np.isnan(v) else f"{'NaN':>20}")
        print("".join(parts))

    print("\n# Stats per station (mean, peak, valid days)\n")
    print(f"{'Station':<32} {'Elev':>5}  {'Valid':>5}  {'Mean':>6}  {'Peak':>6}  {'PeakDay':>11}")
    stats = {}
    for c in targets:
        vals_pairs = [(r["date"], r[c]) for r in rows if not np.isnan(r[c])]
        if not vals_pairs:
            print(f"{STATIONS[c][0]:<32} {STATIONS[c][1]:>5}  {0:>5}  {'-':>6}  {'-':>6}  {'-':>11}")
            stats[c] = {"valid": 0}
            continue
        mean_v = float(np.mean([v for _, v in vals_pairs]))
        peak = max(vals_pairs, key=lambda x: x[1])
        print(
            f"{STATIONS[c][0]:<32} {STATIONS[c][1]:>5}  {len(vals_pairs):>5}  "
            f"{mean_v:>6.2f}  {peak[1]:>6.2f}  {peak[0]:>11}"
        )
        stats[c] = {"valid": len(vals_pairs), "mean_m3s": mean_v, "peak_m3s": peak[1], "peak_date": peak[0]}

    # Mass-balance check: sum upstream stations vs Santa Juana.
    upstream_codes = [
        "03815001", "03815002", "03814003", "03806001", "03804002", "03802001",
    ]
    upstream_total_mean = sum(stats[c].get("mean_m3s", 0.0) for c in upstream_codes if stats[c].get("valid", 0) > 0)
    sj_mean = stats["03820003"].get("mean_m3s", float("nan"))
    print(
        f"\nUpstream tributary sum (mean Q over window): {upstream_total_mean:.2f} m³/s "
        f"({sum(1 for c in upstream_codes if stats[c].get('valid', 0) > 0)} stations with data)"
    )
    print(f"Santa Juana mean Q over window:              {sj_mean:.2f} m³/s")
    print(f"Δ (Santa Juana − upstream sum):              {sj_mean - upstream_total_mean:+.2f} m³/s")
    print(
        "→ The deficit indicates the 2017 event came from sub-basins between\n"
        "  the Tránsito/Carmen confluence and Santa Juana that have NO DGA gauges.\n"
        "  A classic upstream-input/downstream-validation split is NOT POSSIBLE\n"
        "  for this event with these stations."
    )

    # Save summary as JSON for traceability.
    out_dir = Path(__file__).parent / "output"
    out_dir.mkdir(exist_ok=True)
    summary = {
        "window": {"start": WIN_START, "end": WIN_END},
        "stations": {c: {"name": STATIONS[c][0], "elev_m": STATIONS[c][1], **stats[c]} for c in targets},
        "mass_balance": {
            "santa_juana_mean_m3s": sj_mean,
            "upstream_sum_mean_m3s": upstream_total_mean,
            "deficit_m3s": sj_mean - upstream_total_mean,
            "interpretation": "Atacama 2017 was a localised event from sub-basins without DGA gauges.",
        },
    }
    out_path = out_dir / "atacama_2017_basin_stations.json"
    out_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"\nWrote summary: {out_path}")


if __name__ == "__main__":
    main()
