#!/usr/bin/env python3
"""claw-yield-units: convert grain yields between bu/ac, lb/ac, kg/ha, t/ha and g/m2,
adjust to a standard moisture, and scale plot weights to per-area yield, for one
value or a CSV of plots. Every result states its bushel weight and moisture basis.

Usage:
    python scripts/yield_units.py --demo
    python scripts/yield_units.py --value 60 --from bu/ac --crop soybean
    python scripts/yield_units.py --value 11.2 --from t/ha --crop corn --moisture 21 --to-moisture 15.5
    python scripts/yield_units.py --plot-weight 9.8 --weight-unit lb --rows 2 --row-spacing 30 \
        --row-spacing-unit in --length 17.5 --length-unit ft --moisture 14.1 --crop soybean
    python scripts/yield_units.py --input plots.csv --crop corn --weight-unit lb --rows 2 \
        --row-spacing 30 --row-spacing-unit in --length 17.5 --length-unit ft

Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT license.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import platform
import shlex
import sys
from datetime import datetime
from pathlib import Path

SKILL = "claw-yield-units"
VERSION = "1.0.0"
SKILL_DIR = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Helpers shared by all AgriClaw Python scripts (copy unchanged)
# ---------------------------------------------------------------------------


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def save_table(rows: list[dict], out_dir: Path, name: str) -> Path:
    path = out_dir / "tables" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def write_repro_bundle(out_dir: Path, argv: list[str], inputs: list[str] = (),
                       packages: list[str] = ()) -> None:
    """commands.sh, environment.txt and SHA-256 checksums of inputs + outputs."""
    rd = out_dir / "reproducibility"
    rd.mkdir(parents=True, exist_ok=True)
    cmd = " ".join(["python", shlex.quote(str(Path(__file__).resolve()))] +
                   [shlex.quote(a) for a in argv])
    stamp = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
    write_text(rd / "commands.sh",
               f"#!/usr/bin/env bash\n# {SKILL} {VERSION}, run {stamp}\n{cmd}\n")
    env = [f"python {platform.python_version()} ({platform.platform()})"]
    for pkg in packages:
        try:
            env.append(f"{pkg} {importlib.metadata.version(pkg)}")
        except importlib.metadata.PackageNotFoundError:
            env.append(f"{pkg} (not installed)")
    write_text(rd / "environment.txt", "\n".join(env) + "\n")
    files = [Path(p) for p in inputs if p and Path(p).is_file()]
    files += sorted(p for p in out_dir.rglob("*") if p.is_file() and rd not in p.parents)
    lines = []
    for p in files:
        try:
            label = p.resolve().relative_to(out_dir.resolve()).as_posix()
        except ValueError:
            label = p.resolve().as_posix()
        lines.append(f"{sha256_file(p)}  {label}")
    write_text(rd / "checksums.sha256", "\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Exact unit definitions
# ---------------------------------------------------------------------------

LB_KG = 0.45359237                  # avoirdupois pound, exact (NIST SP 811, Appendix B footnote 22)
FT_M = 0.3048                       # foot, exact (international foot; US survey foot retired 2023, 85 FR 62698)
IN_M = 0.0254                       # inch, exact
ACRE_M2 = 43560 * FT_M ** 2         # 4046.8564224 m2
HA_M2 = 10000.0
LB_AC_TO_KG_HA = LB_KG / (ACRE_M2 / HA_M2)   # 1.12085116...

# yield units -> kg/ha (bu/ac handled separately because it needs the bushel weight)
TO_KG_HA = {"kg/ha": 1.0, "t/ha": 1000.0, "mg/ha": 1000.0, "g/m2": 10.0, "lb/ac": LB_AC_TO_KG_HA,
            "cwt/ac": 100.0 * LB_AC_TO_KG_HA}
UNIT_ALIASES = {"kg/ha": "kg/ha", "kg ha-1": "kg/ha", "t/ha": "t/ha", "mg/ha": "t/ha", "metric tons/ha": "t/ha",
                "g/m2": "g/m2", "g/m^2": "g/m2", "lb/ac": "lb/ac", "lb/acre": "lb/ac", "lbs/ac": "lb/ac",
                "cwt/ac": "cwt/ac", "bu/ac": "bu/ac", "bu/acre": "bu/ac", "bushels/acre": "bu/ac"}
WEIGHT_TO_KG = {"g": 0.001, "kg": 1.0, "lb": LB_KG}
LENGTH_TO_M = {"m": 1.0, "cm": 0.01, "in": IN_M, "ft": FT_M}
AREA_TO_M2 = {"m2": 1.0, "ft2": FT_M ** 2, "ha": HA_M2, "ac": ACRE_M2}

# Bushel weights: USDA AH-697 (1992) Tables 5-6, except canola and sunflower (UF/IFAS AG442,
# because AH-697 gives ranges). Standard moistures: UF/IFAS AG442 Table 1; corn, soybean and wheat
# also match Arkansas FSA-1078 Table 1. None = sources disagree or are silent: the user must say.
CROPS = {
    "corn": {"label": "corn (shelled)", "lb_per_bu": 56.0, "std_moisture_pct": 15.5},
    "soybean": {"label": "soybean", "lb_per_bu": 60.0, "std_moisture_pct": 13.0},
    "wheat": {"label": "wheat", "lb_per_bu": 60.0, "std_moisture_pct": 13.5},
    "sorghum": {"label": "grain sorghum", "lb_per_bu": 56.0, "std_moisture_pct": None,
                "note": "standard moisture differs between sources (13% UF/IFAS AG442, 14% Arkansas FSA-1078, "
                        "which also lists 55 lb/bu)"},
    "barley": {"label": "barley", "lb_per_bu": 48.0, "std_moisture_pct": None},
    "oats": {"label": "oats", "lb_per_bu": 32.0, "std_moisture_pct": 13.5},
    "rye": {"label": "rye", "lb_per_bu": 56.0, "std_moisture_pct": None},
    "flaxseed": {"label": "flaxseed", "lb_per_bu": 56.0, "std_moisture_pct": None},
    "rice": {"label": "rough rice", "lb_per_bu": 45.0, "std_moisture_pct": None,
             "note": "US rice is usually reported in cwt (100 lb); Arkansas FSA-1078 lists 13% base moisture"},
    "canola": {"label": "canola", "lb_per_bu": 50.0, "std_moisture_pct": 8.5,
               "note": "AH-697 lists rapeseed at 50-60 lb/bu; 50 lb/bu is from UF/IFAS AG442"},
    "sunflower": {"label": "sunflower (unshelled)", "lb_per_bu": 25.0, "std_moisture_pct": 10.0,
                  "note": "AH-697 lists 24-32 lb/bu; 25 lb/bu is from UF/IFAS AG442. Sunflower is usually "
                          "traded in cwt or lb"},
}


def ascii_safe(text: str) -> str:
    for a, b in (("±", "+/-"), ("²", "2"), ("×", "x"), ("≥", ">="), ("≤", "<=")):
        text = text.replace(a, b)
    return text.encode("ascii", "replace").decode("ascii")


# ---------------------------------------------------------------------------
# Conversions
# ---------------------------------------------------------------------------


def norm_unit(unit: str) -> str:
    u = unit.strip().lower().replace(" ", "")
    if u not in UNIT_ALIASES:
        raise ValueError(f"Unknown yield unit '{unit}'. Use one of: kg/ha, t/ha, g/m2, lb/ac, cwt/ac, bu/ac.")
    return UNIT_ALIASES[u]


def to_kg_ha(value: float, unit: str, lb_per_bu: float | None = None) -> float:
    u = norm_unit(unit)
    if u == "bu/ac":
        if not lb_per_bu:
            raise ValueError("bu/ac needs a bushel weight: give --crop or --lb-per-bu.")
        return value * lb_per_bu * LB_AC_TO_KG_HA
    return value * TO_KG_HA[u]


def from_kg_ha(value_kg_ha: float, unit: str, lb_per_bu: float | None = None) -> float:
    u = norm_unit(unit)
    if u == "bu/ac":
        if not lb_per_bu:
            raise ValueError("bu/ac needs a bushel weight: give --crop or --lb-per-bu.")
        return value_kg_ha / (lb_per_bu * LB_AC_TO_KG_HA)
    return value_kg_ha / TO_KG_HA[u]


def convert(value: float, from_unit: str, to_unit: str, lb_per_bu: float | None = None) -> float:
    return from_kg_ha(to_kg_ha(value, from_unit, lb_per_bu), to_unit, lb_per_bu)


def adjust_moisture(weight: float, moisture_pct: float, target_pct: float) -> float:
    """Weight at a target moisture with the same dry matter (wet basis):
    W_target = W x (100 - M) / (100 - M_target)."""
    for m, name in ((moisture_pct, "moisture"), (target_pct, "target moisture")):
        if not 0 <= m < 100:
            raise ValueError(f"{name} must be a wet-basis percentage between 0 and 100 (got {m}).")
    return weight * (100.0 - moisture_pct) / (100.0 - target_pct)


def plot_area_m2(rows: float | None = None, row_spacing: float | None = None, length: float | None = None,
                 row_spacing_unit: str = "m", length_unit: str = "m", area: float | None = None,
                 area_unit: str = "m2") -> float:
    if area is not None:
        if area_unit not in AREA_TO_M2:
            raise ValueError(f"Unknown area unit '{area_unit}'; use {', '.join(AREA_TO_M2)}.")
        a = area * AREA_TO_M2[area_unit]
    else:
        if None in (rows, row_spacing, length):
            raise ValueError("Give the harvested area (--plot-area) or rows, row spacing and harvested length.")
        for u, table in ((row_spacing_unit, LENGTH_TO_M), (length_unit, LENGTH_TO_M)):
            if u not in table:
                raise ValueError(f"Unknown length unit '{u}'; use {', '.join(LENGTH_TO_M)}.")
        a = rows * row_spacing * LENGTH_TO_M[row_spacing_unit] * length * LENGTH_TO_M[length_unit]
    if a <= 0:
        raise ValueError("Harvested area must be positive.")
    return a


def plot_yield_kg_ha(weight: float, weight_unit: str, area_m2: float) -> float:
    if weight_unit not in WEIGHT_TO_KG:
        raise ValueError(f"Unknown weight unit '{weight_unit}'; use g, kg or lb.")
    if weight < 0:
        raise ValueError("Plot weight cannot be negative.")
    return weight * WEIGHT_TO_KG[weight_unit] / area_m2 * HA_M2


def all_units(kg_ha: float, lb_per_bu: float | None) -> dict:
    out = {"yield_kg_ha": kg_ha, "yield_t_ha": kg_ha / 1000.0, "yield_g_m2": kg_ha / 10.0,
           "yield_lb_ac": kg_ha / LB_AC_TO_KG_HA}
    out["yield_bu_ac"] = kg_ha / (lb_per_bu * LB_AC_TO_KG_HA) if lb_per_bu else None
    return out


# ---------------------------------------------------------------------------
# Basis resolution and checks
# ---------------------------------------------------------------------------


def resolve_basis(crop: str | None, lb_per_bu: float | None, to_moisture: float | None):
    """Bushel weight and target moisture with their sources."""
    info = CROPS.get(crop) if crop else None
    notes = []
    if lb_per_bu is not None:
        bu, bu_src = lb_per_bu, "user (--lb-per-bu)"
        if info and abs(lb_per_bu - info["lb_per_bu"]) > 1e-9:
            notes.append(f"--lb-per-bu {lb_per_bu:g} differs from the {info['lb_per_bu']:g} lb standard for "
                         f"{info['label']}.")
    elif info:
        bu = info["lb_per_bu"]
        bu_src = ("UF/IFAS AG442" if crop in ("canola", "sunflower") else "USDA AH-697")
    else:
        bu, bu_src = None, ""
    if to_moisture is not None:
        tm, tm_src = to_moisture, "user (--to-moisture)"
    elif info and info["std_moisture_pct"] is not None:
        tm = info["std_moisture_pct"]
        tm_src = "UF/IFAS AG442" + ("; Arkansas FSA-1078" if crop in ("corn", "soybean", "wheat") else "")
    else:
        tm, tm_src = None, ""
    if info and info.get("note"):
        notes.append(f"{info['label']}: {info['note']}.")
    return bu, bu_src, tm, tm_src, notes


def moisture_flags(m: float | None) -> list[str]:
    if m is None:
        return []
    if not 5 <= m <= 40:
        return [f"harvest moisture {m:g}% is outside 5-40% (rule of thumb): check it is wet-basis percent"]
    return []


def geometry_flags(area_m2: float, row_spacing_m: float | None, length_m: float | None) -> list[str]:
    flags = []
    if not 0.5 <= area_m2 <= 10000:
        flags.append(f"harvested area {area_m2:.3g} m² is outside 0.5 m²-1 ha (rule of thumb): check units")
    if row_spacing_m is not None and not 0.1 <= row_spacing_m <= 2.0:
        flags.append(f"row spacing {row_spacing_m:.3g} m is outside 0.1-2 m: check --row-spacing-unit")
    if length_m is not None and not 1.0 <= length_m <= 500:
        flags.append(f"harvested length {length_m:.3g} m is outside 1-500 m: check --length-unit")
    return flags


def make_figure(records: list[dict], label_col: str | None, bu: float | None, basis: str, out: Path) -> list[str]:
    """Dot chart of per-area yields (batch and plot modes); flagged records in red."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    vals = [r["yield_kg_ha"] for r in records]
    labels = [str(r.get(label_col, i + 1)) if label_col else str(i + 1) for i, r in enumerate(records)]
    fig, ax = plt.subplots(figsize=(7, max(2.5, 0.28 * len(records) + 1.2)))
    colors = ["#c0392b" if r["flags"] else "#1f5f8b" for r in records]
    ax.scatter(vals, range(len(vals)), c=colors, s=24, zorder=3)
    ax.set_yticks(range(len(vals)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel(f"Yield (kg/ha, {basis})")
    if bu:
        sec = ax.secondary_xaxis("top", functions=(lambda x: x / (bu * LB_AC_TO_KG_HA),
                                                   lambda x: x * bu * LB_AC_TO_KG_HA))
        sec.set_xlabel(f"Yield (bu/ac, {bu:g} lb/bu)")
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    (out / "figures").mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "figures" / "yields.png", dpi=150)
    plt.close(fig)
    return ["figures/yields.png"]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

REFERENCES = [
    "USDA Economic Research Service (1992). Weights, Measures, and Conversion Factors for Agricultural "
    "Commodities and Their Products. Agricultural Handbook No. 697. Tables 5-6. "
    "https://www.ers.usda.gov/publications/pub-details?pubid=41881",
    "Mulvaney, M.J. & Devkota, P. (2020). Adjusting Crop Yield to a Standard Moisture Content. UF/IFAS "
    "Extension SS-AGR-443/AG442. EDIS 2020(3). https://doi.org/10.32473/edis-ag442-2020",
    "Sadaka, S., Atungulu, G.G. & Olatunde, G. Understanding Grain Shrinkage and Expansion. FSA1078, "
    "University of Arkansas Division of Agriculture. https://www.uaex.uada.edu/publications/pdf/FSA-1078.pdf",
    "Thompson, A. & Taylor, B.N. (2008). Guide for the Use of the International System of Units (SI). NIST "
    "Special Publication 811. https://doi.org/10.6028/NIST.SP.811e2008",
    "NIST & NOAA (2020). Deprecation of the United States (U.S.) Survey Foot. Federal Register 85 FR 62698. "
    "https://www.federalregister.gov/documents/2020/10/05/2020-21902/deprecation-of-the-united-states-us-survey-foot",
]


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="yield_units.py", description=__doc__.splitlines()[0])
    ap.add_argument("--demo", action="store_true", help="12 synthetic soybean plots (demo/plot_yields.csv)")
    ap.add_argument("--crop", choices=sorted(CROPS), help="sets bushel weight and standard moisture")
    ap.add_argument("--lb-per-bu", type=float, help="bushel weight in lb (overrides the crop's)")
    ap.add_argument("--moisture", type=float, help="harvest moisture, %% wet basis (single value / plot modes)")
    ap.add_argument("--to-moisture", type=float, help="report yields at this moisture, %% (default: crop standard)")
    ap.add_argument("--value", type=float, help="a yield value to convert")
    ap.add_argument("--from", dest="from_unit", help="unit of --value or of the 'yield' column: kg/ha, t/ha, g/m2, "
                                                     "lb/ac, cwt/ac, bu/ac")
    ap.add_argument("--plot-weight", type=float, help="harvested grain weight of one plot")
    ap.add_argument("--weight-unit", choices=sorted(WEIGHT_TO_KG), help="unit of plot weights")
    ap.add_argument("--rows", type=float, help="number of harvested rows")
    ap.add_argument("--row-spacing", type=float, help="row spacing")
    ap.add_argument("--row-spacing-unit", choices=sorted(LENGTH_TO_M), help="unit of --row-spacing")
    ap.add_argument("--length", type=float, help="harvested row length (after end trimming)")
    ap.add_argument("--length-unit", choices=sorted(LENGTH_TO_M), help="unit of --length")
    ap.add_argument("--plot-area", type=float, help="harvested area (instead of rows x spacing x length)")
    ap.add_argument("--area-unit", choices=sorted(AREA_TO_M2), default="m2", help="unit of --plot-area")
    ap.add_argument("--input", help="CSV: one row per plot (plot_weight or yield column; optional "
                                    "harvest_moisture_pct, rows, row_spacing, length, plot_area)")
    ap.add_argument("--out", default=f"agriclaw_output/{SKILL}", help="output folder")
    return ap


def num_or_none(v):
    v = (v or "").strip() if isinstance(v, str) else v
    if v in ("", None) or (isinstance(v, str) and v.lower() in ("na", "nan")):
        return None
    return float(v)


def process_record(rec: dict, args, bu: float | None, tm: float | None) -> dict:
    """One plot or value -> yields in all units at the target moisture (if any)."""
    flags = []
    m = rec.get("moisture")
    flags += moisture_flags(m)
    area = None
    if rec.get("plot_weight") is not None:
        if not args.weight_unit:
            raise ValueError("--weight-unit (g, kg or lb) is required for plot weights.")
        rs = rec.get("row_spacing")
        ln = rec.get("length")
        if rec.get("plot_area") is None and (rs is not None or ln is not None):
            if not args.row_spacing_unit or not args.length_unit:
                raise ValueError("--row-spacing-unit and --length-unit are required (units are never guessed).")
        area = plot_area_m2(rec.get("rows"), rs, ln, args.row_spacing_unit or "m", args.length_unit or "m",
                            rec.get("plot_area"), args.area_unit)
        flags += geometry_flags(area, rs * LENGTH_TO_M[args.row_spacing_unit] if rs is not None else None,
                                ln * LENGTH_TO_M[args.length_unit] if ln is not None else None)
        kg_ha_raw = plot_yield_kg_ha(rec["plot_weight"], args.weight_unit, area)
    else:
        if not args.from_unit:
            raise ValueError("--from is required to say the unit of the yield value(s).")
        kg_ha_raw = to_kg_ha(rec["value"], args.from_unit, bu)
    if m is not None and tm is not None:
        kg_ha = adjust_moisture(kg_ha_raw, m, tm)
        basis = f"{tm:g}% moisture"
    elif m is not None:
        kg_ha = kg_ha_raw
        basis = f"{m:g}% moisture (as harvested; no target moisture)"
    else:
        kg_ha = kg_ha_raw
        basis = "same as input (moisture not given)"
        if rec.get("value") is not None and args.from_unit and norm_unit(args.from_unit) == "bu/ac" and tm is not None:
            basis += f"; bu/ac are normally quoted at the {tm:g}% standard"
    out = {"area_m2": round(area, 4) if area is not None else "", "harvest_moisture_pct": "" if m is None else m,
           "moisture_basis": basis, "lb_per_bu": "" if bu is None else bu}
    for k, v in all_units(kg_ha, bu).items():
        out[k] = "" if v is None else round(v, 3)
    out["dry_matter_kg_ha"] = round(kg_ha_raw * (100 - m) / 100, 3) if m is not None else ""
    out["flags"] = "; ".join(flags)
    return out


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = build_parser()
    args = ap.parse_args(argv)
    warnings: list[str] = []
    out = Path(args.out)
    try:
        if args.demo:
            args.input = args.input or str(SKILL_DIR / "demo" / "plot_yields.csv")
            args.crop = args.crop or "soybean"
            args.weight_unit = args.weight_unit or "lb"
            args.rows = args.rows or 2
            args.row_spacing = args.row_spacing or 30
            args.row_spacing_unit = args.row_spacing_unit or "in"
            args.length = args.length or 17.5
            args.length_unit = args.length_unit or "ft"
        modes = [x for x in (args.value is not None, args.plot_weight is not None, args.input is not None) if x]
        if len(modes) != 1:
            ap.error("give exactly one of --value, --plot-weight or --input (or --demo)")
        bu, bu_src, tm, tm_src, notes = resolve_basis(args.crop, args.lb_per_bu, args.to_moisture)
        warnings += notes
        if args.crop is None and args.lb_per_bu is None:
            warnings.append("No crop or bushel weight given: bu/ac is not reported.")
        if args.crop in CROPS and CROPS[args.crop]["std_moisture_pct"] is None and args.to_moisture is None:
            if args.moisture is not None or args.input:
                raise ValueError(f"No agreed standard moisture for {CROPS[args.crop]['label']}: give --to-moisture "
                                 "(ask the user which market or trial-network standard applies).")
        records, inputs = [], []
        if args.input:
            p = Path(args.input)
            if not p.is_file():
                raise ValueError(f"Input file not found: {args.input}")
            inputs.append(args.input)
            with p.open(newline="", encoding="utf-8-sig") as f:
                rows = list(csv.DictReader(f))
            if not rows:
                raise ValueError("The input CSV has no rows.")
            cols = set(rows[0])
            if "plot_weight" not in cols and "yield" not in cols:
                raise ValueError(f"The CSV needs a 'plot_weight' or a 'yield' column; found {', '.join(sorted(cols))}.")
            for i, r in enumerate(rows, start=2):
                try:
                    rec = {"moisture": num_or_none(r.get("harvest_moisture_pct")) if "harvest_moisture_pct" in cols
                           else args.moisture}
                    if "plot_weight" in cols:
                        rec["plot_weight"] = num_or_none(r["plot_weight"])
                        if rec["plot_weight"] is None:
                            raise ValueError("plot_weight is missing")
                        rec["rows"] = num_or_none(r.get("rows")) if "rows" in cols else args.rows
                        rec["row_spacing"] = num_or_none(r.get("row_spacing")) if "row_spacing" in cols else args.row_spacing
                        rec["length"] = num_or_none(r.get("length")) if "length" in cols else args.length
                        rec["plot_area"] = num_or_none(r.get("plot_area")) if "plot_area" in cols else args.plot_area
                    else:
                        rec["value"] = num_or_none(r["yield"])
                        if rec["value"] is None:
                            raise ValueError("yield is missing")
                    res = process_record(rec, args, bu, tm)
                except ValueError as e:
                    raise ValueError(f"Row {i}: {e}") from None
                records.append(dict(r, **res))
            if "harvest_moisture_pct" not in cols and args.moisture is None:
                warnings.append("No harvest moisture given: yields are on an unknown moisture basis and were not "
                                "adjusted. Ask for moisture before comparing with other data.")
        else:
            rec = {"moisture": args.moisture}
            if args.value is not None:
                rec["value"] = args.value
            else:
                rec.update(plot_weight=args.plot_weight, rows=args.rows, row_spacing=args.row_spacing,
                           length=args.length, plot_area=args.plot_area)
            records.append(dict({"input": f"{args.value if args.value is not None else args.plot_weight}"}, **process_record(rec, args, bu, tm)))
            if args.moisture is None:
                warnings.append("No harvest moisture given: the result keeps the input's moisture basis (unknown "
                                "unless the user states it).")
        if any(r["flags"] for r in records):
            warnings.append(f"{sum(1 for r in records if r['flags'])} record(s) have plausibility flags (see the "
                            "'flags' column).")
        if bu is not None and tm is not None and abs(tm - (CROPS.get(args.crop or "", {}).get("std_moisture_pct") or tm)) > 1e-9:
            warnings.append(f"Yields are reported at {tm:g}% moisture, not the market standard for this crop; bushels "
                            "are still counted at the standard bushel weight.")

        out.mkdir(parents=True, exist_ok=True)
        save_table(records, out, "yields.csv")
        figures = []
        if args.input or args.plot_weight is not None:
            label_col = next((c for c in ("plot", "plot_id", "id", "entry") if records and c in records[0]), None)
            figures = make_figure(records, label_col, bu, f"{tm:g}% moisture" if tm is not None and
                                  records[0]["harvest_moisture_pct"] != "" else "moisture basis per row", out)
        basis_lines = [f"- Bushel weight: {bu:g} lb/bu ({bu_src})" if bu is not None else
                       "- Bushel weight: none (bu/ac not reported)",
                       f"- Moisture basis: {tm:g}% ({tm_src})" if tm is not None else
                       "- Moisture basis: no target moisture; see each row's `moisture_basis`",
                       f"- 1 lb/ac = {LB_AC_TO_KG_HA:.8f} kg/ha (1 lb = 0.45359237 kg; 1 acre = 4046.8564224 m²)"]
        L = [f"# {SKILL} report", ""]
        if warnings:
            L += ["## Warnings", ""] + [f"- {w}" for w in warnings] + [""]
        L += ["## Basis (always quote this with the numbers)", ""] + basis_lines + [
              f"- Crop: {CROPS[args.crop]['label'] if args.crop else 'not given'}", "",
              "## Results", "",
              "| # | kg/ha | t/ha | lb/ac | bu/ac | Moisture basis | Flags |", "|---|---|---|---|---|---|---|"]
        for i, r in enumerate(records[:50], start=1):
            L.append(f"| {i} | {r['yield_kg_ha']:.1f} | {r['yield_t_ha']:.3f} | {r['yield_lb_ac']:.1f} | "
                     f"{r['yield_bu_ac'] if r['yield_bu_ac'] == '' else format(r['yield_bu_ac'], '.2f')} | "
                     f"{r['moisture_basis']} | {r['flags']} |")
        if len(records) > 50:
            L.append(f"| ... | {len(records) - 50} more rows in tables/yields.csv | | | | | |")
        L += [""] + [f"![{f}]({f})" for f in figures] + ["", "## Method", "",
              "- Moisture adjustment keeps dry matter constant: W_std = W × (100 − M_harvest) / (100 − M_std), with "
              "moisture in % wet basis. Subtracting water and adding back the standard moisture is a common error "
              "(UF/IFAS AG442).",
              "- Plot yield = plot weight / (harvested rows × row spacing × harvested length); use the harvested "
              "length after end trimming, not the planted plot length.",
              "- Bushels are volume units defined by a standard weight at standard moisture; convert to the standard "
              "moisture before converting to bushels.", "",
              "## References", ""] + [f"- {r}" for r in REFERENCES] + ["", f"*{SKILL} {VERSION}*"]
        write_text(out / "report.md", "\n".join(L) + "\n")
        write_repro_bundle(out, argv, inputs=inputs, packages=["matplotlib"])
    except ValueError as e:
        print(ascii_safe(f"Error: {e}"), file=sys.stderr)
        return 1
    r0 = records[0]
    print(ascii_safe(f"Done: {len(records)} record(s). First: {r0['yield_kg_ha']:.1f} kg/ha, {r0['yield_lb_ac']:.1f} lb/ac"
                     + (f", {r0['yield_bu_ac']:.2f} bu/ac at {bu:g} lb/bu" if r0["yield_bu_ac"] != "" else "")
                     + f" ({r0['moisture_basis']}). Report: {(out / 'report.md').resolve()}"))
    for w in warnings:
        print(ascii_safe(f"Warning: {w}"), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
