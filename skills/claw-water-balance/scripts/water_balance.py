#!/usr/bin/env python3
"""claw-water-balance: FAO-56 single-crop-coefficient daily root-zone soil water
balance: Kc curve, ETc, readily and total available water, water-stress coefficient
Ks, adjusted ETc, deep percolation, stress days and optional irrigation scheduling.

Usage:
    python scripts/water_balance.py --demo
    python scripts/water_balance.py --weather daily_weather.csv --crop maize_grain \
        --plant-date 2025-04-25 --soil-dir results/soil --root-depth-m 1.2
    python scripts/water_balance.py --weather w.csv --crop soybean --plant-date 2025-05-20 \
        --theta-fc 0.32 --theta-wp 0.15 --irrigate

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
from datetime import date, datetime, timedelta
from pathlib import Path

SKILL = "claw-water-balance"
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


def ascii_safe(text: str) -> str:
    """Console-safe text (Windows consoles may not encode symbols such as ° or ≥)."""
    for a, b in (("°", "deg"), ("≥", ">="), ("≤", "<="), ("θ", "theta"), ("−", "-"), ("–", "-"), ("×", "x"),
                 ("²", "2"), ("±", "+/-")):
        text = text.replace(a, b)
    return text.encode("ascii", "replace").decode("ascii")


# ---------------------------------------------------------------------------
# FAO-56 crop parameters (Allen et al. 1998): Table 12 (Kc, max height h),
# Table 22 (max effective root depth range, p for ETc ~ 5 mm/day), Table 11
# (one row of stage lengths: initial, development, mid-season, late season, days).
# Where FAO gives a range without guidance, the midpoint is used (noted).
# ---------------------------------------------------------------------------

CROPS = {
    "maize_grain": dict(name="Maize, field (grain)", kc_ini=0.30, kc_mid=1.20, kc_end=0.35, h=2.0,
                        zr=(1.0, 1.7), p=0.55, stages=(30, 40, 50, 30),
                        stage_row="Table 11: April planting, Spain (spring/summer) and California",
                        note="Kc_end 0.35 is for grain dried in the field (~18% moisture); use 0.60 for "
                             "harvest at high grain moisture (Table 12 footnote 11)."),
    "soybean": dict(name="Soybeans", kc_ini=0.40, kc_mid=1.15, kc_end=0.50, h=0.75,
                    zr=(0.6, 1.3), p=0.50, stages=(20, 35, 60, 25),
                    stage_row="Table 11: May planting, Central USA (development listed as 30/35 days)",
                    note="h = 0.75 m is the midpoint of Table 12's 0.5-1.0 m."),
    "spring_wheat": dict(name="Spring wheat", kc_ini=0.30, kc_mid=1.15, kc_end=0.25, h=1.0,
                         zr=(1.0, 1.5), p=0.55, stages=(20, 25, 60, 30),
                         stage_row="Table 11: barley/oats/wheat, March/April planting, 35-45° latitude",
                         note="Kc_end 0.25 for machine harvest; 0.4 for hand harvest (Table 12 footnote 10)."),
    "winter_wheat": dict(name="Winter wheat", kc_ini=0.40, kc_mid=1.15, kc_end=0.25, h=1.0,
                         zr=(1.5, 1.8), p=0.55, stages=(160, 75, 75, 25),
                         stage_row="Table 11: October planting, Idaho, USA",
                         note="Kc_ini 0.4 is for frozen soils; use 0.7 where soils do not freeze (Table 12)."),
    "barley": dict(name="Barley", kc_ini=0.30, kc_mid=1.15, kc_end=0.25, h=1.0,
                   zr=(1.0, 1.5), p=0.55, stages=(20, 25, 60, 30),
                   stage_row="Table 11: barley/oats/wheat, March/April planting, 35-45° latitude", note=""),
    "sorghum_grain": dict(name="Sorghum, grain", kc_ini=0.30, kc_mid=1.05, kc_end=0.55, h=1.5,
                          zr=(1.0, 2.0), p=0.55, stages=(20, 35, 40, 30),
                          stage_row="Table 11: May/June planting, USA, Pakistan, Mediterranean (listed total 130 d)",
                          note="Kc_mid 1.05 and h 1.5 m are midpoints of 1.00-1.10 and 1-2 m."),
    "cotton": dict(name="Cotton", kc_ini=0.35, kc_mid=1.175, kc_end=0.60, h=1.35,
                   zr=(1.0, 1.7), p=0.65, stages=(30, 50, 55, 45),
                   stage_row="Table 11: April planting, Texas",
                   note="Kc_mid, Kc_end and h are midpoints of 1.15-1.20, 0.70-0.50 and 1.2-1.5 m."),
    "dry_bean": dict(name="Beans, dry", kc_ini=0.40, kc_mid=1.15, kc_end=0.35, h=0.4,
                     zr=(0.6, 0.9), p=0.45, stages=(20, 30, 40, 20),
                     stage_row="Table 11: May/June planting, continental climates", note=""),
    "sunflower": dict(name="Sunflower", kc_ini=0.35, kc_mid=1.075, kc_end=0.35, h=2.0,
                      zr=(0.8, 1.5), p=0.45, stages=(25, 35, 45, 25),
                      stage_row="Table 11: April/May planting, Mediterranean and California",
                      note="Kc_mid 1.075 is the midpoint of 1.0-1.15; lower values are for rainfed, less "
                           "dense stands (Table 12 footnote 9)."),
    "potato": dict(name="Potato", kc_ini=0.50, kc_mid=1.15, kc_end=0.75, h=0.6,
                   zr=(0.4, 0.6), p=0.35, stages=(25, 30, 45, 30),
                   stage_row="Table 11: May planting, continental climate",
                   note="Kc_end is about 0.40 for long-season potatoes with vine kill (Table 12 footnote 4)."),
}
STAGE_NAMES = ("initial", "development", "mid-season", "late season")

DEMO = {"weather": "daily_weather_champaign_2025.csv", "soil_dir": "soil_drummer", "crop": "maize_grain",
        "plant_date": "2025-04-25", "root_depth_m": 1.2}


# ---------------------------------------------------------------------------
# FAO-56 equations
# ---------------------------------------------------------------------------


def taw_from_theta(theta_fc: float, theta_wp: float, zr_m: float) -> float:
    """FAO-56 eq. 82: TAW = 1000 (theta_FC - theta_WP) Zr, mm."""
    if not 0 < theta_wp < theta_fc < 0.8:
        raise ValueError("Need 0 < theta_wp < theta_fc < 0.8 (volumetric, m3/m3); e.g. 0.32 and 0.15. If you "
                         "have percent values, divide by 100.")
    return 1000.0 * (theta_fc - theta_wp) * zr_m


def kc_climate_adjust(kc_tab: float, u2: float, rh_min: float, h: float) -> float:
    """FAO-56 eq. 62 (Kc_mid) and eq. 65 (Kc_end), with the stated limits on u2, RHmin and h."""
    u2 = min(max(u2, 1.0), 6.0)
    rh_min = min(max(rh_min, 20.0), 80.0)
    h = min(max(h, 0.1), 10.0)
    return kc_tab + (0.04 * (u2 - 2.0) - 0.004 * (rh_min - 45.0)) * (h / 3.0) ** 0.3


def kc_on_day(i: int, stages, kc_ini: float, kc_mid: float, kc_end: float):
    """FAO-56 eq. 66: Kc on day i (1 = planting day). Returns (Kc, stage) or (None, '') after the season."""
    l_ini, l_dev, l_mid, l_late = stages
    if i <= l_ini:
        return kc_ini, STAGE_NAMES[0]
    if i <= l_ini + l_dev:
        return kc_ini + (i - l_ini) / l_dev * (kc_mid - kc_ini), STAGE_NAMES[1]
    if i <= l_ini + l_dev + l_mid:
        return kc_mid, STAGE_NAMES[2]
    if i <= l_ini + l_dev + l_mid + l_late:
        return kc_mid + (i - (l_ini + l_dev + l_mid)) / l_late * (kc_end - kc_mid), STAGE_NAMES[3]
    return None, ""


def p_adjusted(p_table: float, etc_mm: float) -> float:
    """FAO-56 Table 22 footnote 2: p = p_table + 0.04 (5 - ETc), limited to 0.1-0.8."""
    return min(max(p_table + 0.04 * (5.0 - etc_mm), 0.1), 0.8)


def ks_from_depletion(dr: float, taw: float, raw: float) -> float:
    """FAO-56 eq. 84: Ks = (TAW - Dr) / (TAW - RAW) when Dr > RAW, else 1."""
    if dr <= raw:
        return 1.0
    return max(0.0, min(1.0, (taw - dr) / (taw - raw)))


def run_balance(days: list[dict], taw: float, p_table: float, dr0: float = 0.0, p_adjust: bool = True,
                auto_irrigate: bool = False, mad: float | None = None,
                min_rain_frac_et0: float = 0.0) -> list[dict]:
    """Daily root-zone depletion (FAO-56 eqs. 83-86) for days with keys
    kc, et0_mm, precip_mm and optional runoff_mm, irrigation_mm, capillary_rise_mm.

    Timing (FAO-56 Example 38): rain, irrigation and capillary rise enter early in the
    day; water above field capacity drains the same day (DP); Ks is then computed from
    the depletion after wetting and the day's ETc_adj is removed. ETc_adj is limited
    so that depletion never exceeds TAW, which keeps the water balance closed.
    """
    if taw <= 0:
        raise ValueError("TAW must be positive.")
    if not 0 <= dr0 <= taw:
        raise ValueError(f"Initial depletion {dr0} mm must be within 0..TAW ({taw:.1f} mm).")
    out, dr_prev = [], dr0
    for d in days:
        etc = d["kc"] * d["et0_mm"]
        p = p_adjusted(p_table, etc) if p_adjust else p_table
        raw = p * taw
        precip = d["precip_mm"]
        runoff = d.get("runoff_mm") or 0.0
        if runoff > precip + 1e-9:
            raise ValueError(f"{d['date']}: runoff ({runoff} mm) exceeds precipitation ({precip} mm).")
        p_net = precip - runoff
        ignored = p_net if (min_rain_frac_et0 > 0 and p_net < min_rain_frac_et0 * d["et0_mm"]) else 0.0
        p_eff = p_net - ignored
        cr = d.get("capillary_rise_mm") or 0.0
        irr = d.get("irrigation_mm") or 0.0
        auto = 0.0
        threshold = (mad if mad is not None else p) * taw
        if auto_irrigate and dr_prev > 0 and dr_prev >= threshold - 1e-9:
            auto = dr_prev           # net depth that refills the root zone to field capacity
            irr += auto
        water_in = p_eff + irr + cr
        dp = max(0.0, water_in - dr_prev)
        dr_start = dr_prev - water_in + dp
        ks = ks_from_depletion(dr_start, taw, raw)
        etc_adj = ks * etc
        limited = etc_adj > taw - dr_start + 1e-12
        if limited:
            etc_adj = taw - dr_start
        dr = dr_start + etc_adj
        out.append({
            "date": d["date"], "kc": d["kc"], "et0_mm": d["et0_mm"], "etc_mm": etc, "p": p, "taw_mm": taw,
            "raw_mm": raw, "precip_mm": precip, "runoff_mm": runoff, "precip_ignored_mm": ignored,
            "eff_precip_mm": p_eff, "irrigation_mm": irr, "auto_irrigation_mm": auto, "capillary_rise_mm": cr,
            "deep_percolation_mm": dp, "dr_start_mm": dr_prev, "dr_after_wetting_mm": dr_start, "ks": ks,
            "etc_adj_mm": etc_adj, "dr_end_mm": dr, "stressed": ks < 1.0, "et_limited_by_storage": limited})
        dr_prev = dr
    return out


def mass_balance_residual(rows: list[dict], dr0: float) -> float:
    """(inputs - ETc_adj - DP) - (Dr0 - Dr_end); zero when the balance closes."""
    inputs = sum(r["eff_precip_mm"] + r["irrigation_mm"] + r["capillary_rise_mm"] for r in rows)
    outputs = sum(r["etc_adj_mm"] + r["deep_percolation_mm"] for r in rows)
    return (inputs - outputs) - (dr0 - rows[-1]["dr_end_mm"])


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

MISSING = {"", "na", "nan", "null", "none"}
WEATHER_REQUIRED = ("date", "precip_mm", "et0_mm")
WEATHER_OPTIONAL = ("rh_min_pct", "wind_2m_m_s", "runoff_mm", "irrigation_mm", "capillary_rise_mm")


def read_weather(path: str) -> list[dict]:
    p = Path(path)
    if not p.is_file():
        raise ValueError(f"Weather file not found: {path}")
    with p.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        cols = [c.strip() for c in (reader.fieldnames or [])]
        raw = list(reader)
    missing = [c for c in WEATHER_REQUIRED if c not in cols]
    if missing:
        raise ValueError(f"The weather file needs columns {', '.join(WEATHER_REQUIRED)} (mm per day); missing "
                         f"{', '.join(missing)}. claw-weather's tables/daily_weather.csv has them. Found: {', '.join(cols)}.")
    rows = []
    for i, r in enumerate(raw, start=2):
        r = {k.strip(): v for k, v in r.items()}
        try:
            d = date.fromisoformat(r["date"].strip())
        except ValueError:
            raise ValueError(f"Row {i}: date '{r['date']}' is not YYYY-MM-DD.") from None
        row = {"date": d}
        for c in WEATHER_REQUIRED[1:] + WEATHER_OPTIONAL:
            v = (r.get(c) or "").strip()
            if v.lower() in MISSING:
                row[c] = None
            else:
                try:
                    row[c] = float(v)
                except ValueError:
                    raise ValueError(f"Row {i}: {c} '{v}' is not numeric.") from None
        rows.append(row)
    rows.sort(key=lambda r: r["date"])
    for c, lo, hi in (("precip_mm", 0, 1000), ("et0_mm", 0, 20), ("rh_min_pct", 0, 100), ("wind_2m_m_s", 0, 60)):
        bad = [r for r in rows if r.get(c) is not None and not lo <= r[c] <= hi]
        if bad:
            raise ValueError(f"{c} outside {lo}..{hi} on {bad[0]['date']} ({bad[0][c]}): check units "
                             "(mm/day, %, m/s at 2 m) and missing-value codes.")
    return rows


def read_soil_dir(soil_dir: str) -> tuple[list[tuple[float, float]], dict]:
    """claw-soil output: tables/aws_profile.csv (depth_cm, aws_mm) and tables/mapunit_summary.csv."""
    d = Path(soil_dir)
    prof_path = d / "tables" / "aws_profile.csv"
    if not prof_path.is_file():
        prof_path = d / "aws_profile.csv"
    if not prof_path.is_file():
        raise ValueError(f"No aws_profile.csv in {soil_dir} (expected the output folder of claw-soil).")
    with prof_path.open(newline="", encoding="utf-8-sig") as f:
        prof = [(float(r["depth_cm"]), float(r["aws_mm"]), float(r.get("max_described_depth_cm") or "nan"))
                for r in csv.DictReader(f) if r.get("aws_mm") not in (None, "")]
    if not prof:
        raise ValueError(f"{prof_path} has no available-water values.")
    summary = {}
    for cand in (d / "tables" / "mapunit_summary.csv", d / "mapunit_summary.csv"):
        if cand.is_file():
            with cand.open(newline="", encoding="utf-8-sig") as f:
                summary = next(csv.DictReader(f), {})
            break
    summary["max_described_depth_cm"] = prof[0][2]
    return [(0.0, 0.0)] + [(a, b) for a, b, _ in prof], summary


def interp_aws(profile: list[tuple[float, float]], depth_cm: float) -> float:
    pts = sorted(profile)
    if depth_cm > pts[-1][0]:
        raise ValueError(f"Root depth {depth_cm:g} cm is deeper than the AWS profile ({pts[-1][0]:g} cm).")
    for (d0, a0), (d1, a1) in zip(pts, pts[1:]):
        if d0 <= depth_cm <= d1:
            return a0 + (a1 - a0) * (depth_cm - d0) / (d1 - d0) if d1 > d0 else a0
    return pts[0][1]


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def runs(flags: list[bool]) -> int:
    best = cur = 0
    for f in flags:
        cur = cur + 1 if f else 0
        best = max(best, cur)
    return best


def make_figures(rows: list[dict], stage_dates: list[tuple[str, date]], out: Path) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    (out / "figures").mkdir(parents=True, exist_ok=True)
    dates = [r["date"] for r in rows]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True, gridspec_kw={"height_ratios": [3, 1.2]})
    a1.plot(dates, [r["dr_end_mm"] for r in rows], color="#8b4513", lw=2, label="root-zone depletion Dr")
    a1.plot(dates, [r["raw_mm"] for r in rows], color="#e08a1e", ls="--", lw=1.2, label="RAW (stress threshold)")
    a1.plot(dates, [r["taw_mm"] for r in rows], color="#555555", ls=":", lw=1.2, label="TAW (wilting point)")
    stressed = [r for r in rows if r["stressed"]]
    if stressed:
        a1.scatter([r["date"] for r in stressed], [r["dr_end_mm"] for r in stressed], s=8, color="#c0392b",
                   zorder=3, label="Ks < 1 (water stress)")
    a1.invert_yaxis()
    a1.set_ylabel("Depletion below field capacity (mm)")
    a1.legend(fontsize=7, loc="lower left")
    a1.grid(alpha=0.3)
    for name, d in stage_dates:
        a1.axvline(d, color="#999999", lw=0.6)
        a1.text(d, 0, f" {name}", fontsize=6, va="top", rotation=90, color="#666666")
    a2.bar(dates, [r["eff_precip_mm"] for r in rows], width=1.0, color="#4a90c2", label="effective precipitation")
    a2.bar(dates, [r["irrigation_mm"] for r in rows], width=1.0, color="#2f6b3a",
           bottom=[r["eff_precip_mm"] for r in rows], label="irrigation")
    a2.set_ylabel("mm/day")
    a2.legend(fontsize=7)
    a2.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(out / "figures" / "depletion.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(dates, [r["kc"] for r in rows], color="#2f6b3a", lw=2, label="Kc (eq. 66)")
    ax.plot(dates, [r["ks"] * r["kc"] for r in rows], color="#c0392b", lw=1, label="Ks × Kc")
    ax.set_ylabel("Coefficient (-)")
    ax2 = ax.twinx()
    ax2.plot(dates, [r["et0_mm"] for r in rows], color="#999999", lw=0.8, label="ET0")
    ax2.plot(dates, [r["etc_adj_mm"] for r in rows], color="#1f5f8b", lw=0.8, label="ETc adj")
    ax2.set_ylabel("mm/day")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper left")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax.grid(alpha=0.3)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(out / "figures" / "kc_et.png", dpi=150)
    plt.close(fig)
    return ["figures/depletion.png", "figures/kc_et.png"]


REFERENCES = [
    "Allen, R.G., Pereira, L.S., Raes, D. & Smith, M. (1998). Crop evapotranspiration: guidelines for computing "
    "crop water requirements. FAO Irrigation and Drainage Paper 56. FAO, Rome. ISBN 92-5-104219-5. "
    "https://www.fao.org/4/x0490e/x0490e00.htm (Chapter 6: Tables 11-12, eqs. 62, 65, 66; Chapter 8: eqs. 82-88, "
    "Table 22, Examples 36-38)",
    "Soil Survey Staff, NRCS, USDA. Soil Survey Geographic (SSURGO) Database (available water storage), via "
    "claw-soil. https://sdmdataaccess.sc.egov.usda.gov/",
]


def write_report(out: Path, ctx: dict) -> None:
    s = ctx["summary"]
    L = [f"# {SKILL} report", ""]
    if ctx["warnings"]:
        L += ["## Warnings", ""] + [f"- {w}" for w in ctx["warnings"]] + [""]
    L += ["## Set-up", "",
          f"- Crop: {ctx['crop_name']}; planting {ctx['plant']}; stage lengths {ctx['stages']} days "
          f"(initial, development, mid-season, late) = season end {ctx['season_end']}",
          f"- Weather: `{ctx['weather']}`; days simulated: {s['days']} ({s['first_date']} to {s['last_date']})",
          f"- Root depth Zr = {ctx['zr']:g} m; TAW = {ctx['taw']:.1f} mm ({ctx['taw_source']})",
          f"- p (Table 22) = {ctx['p_table']:g}" + (", adjusted daily for ETc (0.1-0.8)" if ctx['p_adjust'] else ", constant"),
          f"- Initial depletion = {ctx['dr0']:.1f} mm" + (" (field capacity assumed)" if ctx['dr0'] == 0 else ""),
          f"- Irrigation: {ctx['irrigation_desc']}", "",
          "## Crop coefficients", "", "| Parameter | Value | Source |", "|---|---|---|"]
    for p in ctx["params"]:
        L.append(f"| {p['parameter']} | {p['value']} | {p['source']} |")
    L += ["", "## Season water balance", "", "| Term | mm |", "|---|---|",
          f"| Precipitation | {s['precip_mm']:.1f} |",
          f"| Runoff (given) | {s['runoff_mm']:.1f} |",
          f"| Small rain ignored | {s['precip_ignored_mm']:.1f} |",
          f"| Irrigation (net) | {s['irrigation_mm']:.1f} |",
          f"| Capillary rise (given) | {s['capillary_rise_mm']:.1f} |",
          f"| Crop ET without stress, ETc | {s['etc_mm']:.1f} |",
          f"| Crop ET with stress, ETc adj | {s['etc_adj_mm']:.1f} ({100 * s['etc_adj_mm'] / s['etc_mm']:.0f}% of ETc) |",
          f"| Deep percolation | {s['deep_percolation_mm']:.1f} |",
          f"| Depletion at start / end | {ctx['dr0']:.1f} / {s['final_depletion_mm']:.1f} |",
          f"| Mass-balance residual | {s['mass_balance_residual_mm']:.2e} |", "",
          f"**Water-stress days (Ks < 1): {s['stress_days']}** (first {s['first_stress_date'] or 'none'}; longest run "
          f"{s['longest_stress_run_days']} days). Stress days by stage: "
          + ", ".join(f"{k} {v}" for k, v in ctx["stress_by_stage"].items()) + ".", ""]
    if ctx["events"]:
        L += ["## Irrigation events (net depth)", "", "| Date | Net irrigation (mm) | Depletion before (mm) |", "|---|---|---|"]
        L += [f"| {e['date']} | {e['net_irrigation_mm']:.1f} | {e['depletion_before_mm']:.1f} |" for e in ctx["events"][:40]]
        L += ["", "Net depths: divide by the application efficiency of the system to get gross depths.", ""]
    L += ["## Figures", ""] + [f"![{f}]({f})" for f in ctx["figures"]] + [""]
    L += ["## Assumptions and limitations", "",
          "- Single crop coefficient (FAO-56 Chapter 6): soil evaporation is averaged into Kc, so day-to-day ET "
          "after rain is smoothed; Ks reduces Kc as a whole (eq. 81).",
          "- Runoff is zero unless a `runoff_mm` column is given; capillary rise is zero unless `capillary_rise_mm` is "
          "given (FAO-56: valid when the water table is more than ~1 m below the root zone).",
          "- Rain, irrigation and capillary rise enter early in the day and water above field capacity drains the same "
          "day (FAO-56 Example 38). FAO-56 eq. 88 subtracts the day's ET before computing drainage, which gives less "
          "deep percolation on wetting days.",
          "- Root depth is constant from planting; FAO-56 Example 38 shows root growth, which is not modeled.",
          "- Kc values are FAO-56 typical values for standard (non-stressed, well-managed) conditions; local "
          "calibration of stage lengths and Kc is preferred.", ""]
    L += ["## References", ""] + [f"- {r}" for r in REFERENCES] + ["", f"*{SKILL} {VERSION}*"]
    write_text(out / "report.md", "\n".join(L) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="water_balance.py", description=__doc__.splitlines()[0])
    ap.add_argument("--demo", action="store_true", help="rainfed maize near Champaign, IL, 2025, on Drummer soil")
    ap.add_argument("--weather", help="daily CSV with date, precip_mm, et0_mm (e.g. claw-weather tables/daily_weather.csv)")
    ap.add_argument("--crop", choices=sorted(CROPS) + ["custom"], help="FAO-56 crop (custom needs all Kc, stages, Zr, p)")
    ap.add_argument("--plant-date", help="planting (or green-up) date, YYYY-MM-DD; day 1 of the season")
    ap.add_argument("--stage-days", help="initial,development,mid,late stage lengths in days (default FAO-56 Table 11)")
    ap.add_argument("--kc-ini", type=float)
    ap.add_argument("--kc-mid", type=float, help="Table 12 value before climate adjustment")
    ap.add_argument("--kc-end", type=float, help="Table 12 value before climate adjustment")
    ap.add_argument("--crop-height-m", type=float, help="mean maximum crop height (for eqs. 62 and 65)")
    ap.add_argument("--no-climate-adjust", action="store_true", help="skip the RHmin/wind adjustment of Kc_mid/Kc_end")
    ap.add_argument("--root-depth-m", type=float, help="effective root depth Zr (default: Table 22 range end)")
    ap.add_argument("--p", type=float, help="depletion fraction p for ETc ~ 5 mm/day (default Table 22)")
    ap.add_argument("--no-p-adjust", action="store_true", help="keep p constant instead of adjusting it for ETc")
    ap.add_argument("--taw-mm", type=float, help="total available water in the root zone, mm")
    ap.add_argument("--theta-fc", type=float, help="volumetric water content at field capacity, m3/m3")
    ap.add_argument("--theta-wp", type=float, help="volumetric water content at wilting point, m3/m3")
    ap.add_argument("--soil-dir", help="claw-soil output folder (uses tables/aws_profile.csv)")
    ap.add_argument("--initial-depletion-mm", type=float, default=0.0, help="Dr at the start (0 = field capacity)")
    ap.add_argument("--irrigate", action="store_true", help="schedule irrigation when depletion reaches MAD")
    ap.add_argument("--mad", type=float, help="management-allowed depletion, fraction of TAW (default: p)")
    ap.add_argument("--ignore-rain-below-et0-frac", type=float, default=0.0,
                    help="ignore daily rain smaller than this fraction of ET0 (FAO-56 suggests ~0.2; default 0)")
    ap.add_argument("--out", default=f"agriclaw_output/{SKILL}", help="output folder")
    return ap


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = build_parser()
    args = ap.parse_args(argv)
    warnings: list[str] = []
    try:
        if args.demo:
            args.weather = args.weather or str(SKILL_DIR / "demo" / DEMO["weather"])
            args.soil_dir = args.soil_dir or (None if (args.taw_mm or args.theta_fc) else
                                              str(SKILL_DIR / "demo" / DEMO["soil_dir"]))
            args.crop = args.crop or DEMO["crop"]
            args.plant_date = args.plant_date or DEMO["plant_date"]
            args.root_depth_m = args.root_depth_m or DEMO["root_depth_m"]
        for flag, val in (("--weather", args.weather), ("--crop", args.crop), ("--plant-date", args.plant_date)):
            if not val:
                ap.error(f"{flag} is required")
        plant = date.fromisoformat(args.plant_date)
        crop = dict(CROPS.get(args.crop) or {})
        if args.crop == "custom":
            need = [f for f, v in (("--kc-ini", args.kc_ini), ("--kc-mid", args.kc_mid), ("--kc-end", args.kc_end),
                                   ("--stage-days", args.stage_days), ("--root-depth-m", args.root_depth_m),
                                   ("--p", args.p)) if v is None]
            if need:
                ap.error(f"--crop custom needs {', '.join(need)}")
            crop = dict(name="custom crop", h=args.crop_height_m, zr=(args.root_depth_m, args.root_depth_m),
                        stage_row="user", note="")
        src = "FAO-56 Table 12"
        kc_ini = args.kc_ini if args.kc_ini is not None else crop["kc_ini"]
        kc_mid_tab = args.kc_mid if args.kc_mid is not None else crop["kc_mid"]
        kc_end_tab = args.kc_end if args.kc_end is not None else crop["kc_end"]
        h = args.crop_height_m if args.crop_height_m is not None else crop.get("h")
        if args.stage_days:
            try:
                stages = tuple(int(x) for x in args.stage_days.split(","))
            except ValueError:
                ap.error("--stage-days must be four integers, e.g. 30,40,50,30")
            if len(stages) != 4 or min(stages) < 1:
                ap.error("--stage-days must be four positive integers (initial,development,mid,late)")
            stage_src = "user"
        else:
            stages = crop["stages"]
            stage_src = f"FAO-56 {crop['stage_row']}"
            warnings.append(f"Stage lengths {stages} are FAO-56 Table 11 generic values ({crop['stage_row']}). Replace "
                            "them with local observations (--stage-days); FAO-56 shows that mid-season length errors "
                            "change seasonal ET substantially.")
        p_table = args.p if args.p is not None else crop["p"]
        if not 0.05 <= p_table <= 0.9:
            raise ValueError("--p must be a fraction between 0.05 and 0.9.")
        if args.root_depth_m is not None:
            zr, zr_src = args.root_depth_m, "user"
        else:
            zr = crop["zr"][0] if args.irrigate else crop["zr"][1]
            zr_src = (f"FAO-56 Table 22 range {crop['zr'][0]}-{crop['zr'][1]} m: "
                      + ("smaller value, for irrigation scheduling" if args.irrigate
                         else "larger value, for rainfed / stress modeling"))
            warnings.append(f"Root depth {zr} m was taken from FAO-56 Table 22 ({zr_src}); restrictive layers, "
                            "shallow water tables or young crops have shallower roots.")
        if not 0.1 <= zr <= 3.0:
            raise ValueError("--root-depth-m must be between 0.1 and 3 m.")

        # total available water
        sources = [x for x in (args.taw_mm is not None, args.theta_fc is not None or args.theta_wp is not None,
                               args.soil_dir is not None) if x]
        if len(sources) != 1:
            raise ValueError("Give exactly one soil water source: --taw-mm, or --theta-fc with --theta-wp, or "
                             "--soil-dir (claw-soil output). Ask the user; never guess soil water holding.")
        soil_summary = {}
        inputs = [args.weather]
        if args.taw_mm is not None:
            taw, taw_src = args.taw_mm, "user --taw-mm"
        elif args.soil_dir:
            profile, soil_summary = read_soil_dir(args.soil_dir)
            maxd = soil_summary.get("max_described_depth_cm")
            if maxd == maxd and maxd and zr * 100 > maxd:
                warnings.append(f"Root depth {zr * 100:g} cm is below the deepest described horizon ({maxd:g} cm); "
                                "no water is counted below it.")
            taw = interp_aws(profile, zr * 100.0)
            taw_src = (f"SSURGO available water storage to {zr * 100:g} cm via claw-soil "
                       f"({soil_summary.get('muname', 'map unit')})")
            inputs.append(str(Path(args.soil_dir) / "tables" / "aws_profile.csv"))
        else:
            if args.theta_fc is None or args.theta_wp is None:
                raise ValueError("Give both --theta-fc and --theta-wp.")
            taw = taw_from_theta(args.theta_fc, args.theta_wp, zr)
            taw_src = f"FAO-56 eq. 82 with θFC {args.theta_fc}, θWP {args.theta_wp}, Zr {zr} m"
        if taw <= 0:
            raise ValueError("TAW must be positive.")
        wt = soil_summary.get("water_table_min_cm")
        if wt not in (None, ""):
            wt = float(wt)
            if wt < zr * 100 + 100:
                warnings.append(f"SSURGO records a water table as shallow as {wt:g} cm, less than 1 m below the "
                                f"{zr * 100:g} cm root zone. Capillary rise is set to zero, which FAO-56 accepts only "
                                "for deeper water tables; ask whether the field is tile-drained, and treat stress "
                                "estimates as conservative (too high) if it is not.")

        weather = read_weather(args.weather)
        season_end = plant + timedelta(days=sum(stages) - 1)
        if weather[0]["date"] > plant:
            raise ValueError(f"The weather file starts on {weather[0]['date']}, after planting ({plant}).")
        by_date = {w["date"]: w for w in weather}
        last = min(season_end, weather[-1]["date"])
        if last < season_end:
            warnings.append(f"Weather ends on {last}, before the season end ({season_end}); results cover planting to "
                            f"{last} only and nothing is extrapolated.")
        days, missing = [], []
        d = plant
        while d <= last:
            w = by_date.get(d)
            if w is None or w["precip_mm"] is None or w["et0_mm"] is None:
                missing.append(d)
            else:
                days.append(w)
            d += timedelta(days=1)
        if missing:
            raise ValueError(f"{len(missing)} season day(s) lack precipitation or ET0 (first {missing[0]}). The water "
                             "balance cannot skip days; fill them from another source and rerun.")

        # Kc_mid / Kc_end climate adjustment (eqs. 62, 65) from mid- and late-season means
        l_ini, l_dev, l_mid, _ = stages
        mid_days = [w for w in days if l_ini + l_dev < (w["date"] - plant).days + 1 <= l_ini + l_dev + l_mid]
        late_days = [w for w in days if (w["date"] - plant).days + 1 > l_ini + l_dev + l_mid]
        params = [{"parameter": "Kc_ini", "value": f"{kc_ini:.3f}", "source": "user" if args.kc_ini is not None else src}]
        kc_mid, kc_end = kc_mid_tab, kc_end_tab

        def stage_means(ds):
            u = [w["wind_2m_m_s"] for w in ds if w.get("wind_2m_m_s") is not None]
            r = [w["rh_min_pct"] for w in ds if w.get("rh_min_pct") is not None]
            return (sum(u) / len(u) if u else None, sum(r) / len(r) if r else None)

        if args.no_climate_adjust:
            warnings.append("Kc_mid and Kc_end were not adjusted for climate (--no-climate-adjust).")
        elif h is None:
            warnings.append("No crop height: Kc_mid and Kc_end were not adjusted for climate.")
        else:
            u2m, rhm = stage_means(mid_days)
            if u2m is not None and rhm is not None:
                kc_mid = kc_climate_adjust(kc_mid_tab, u2m, rhm, h)
                params.append({"parameter": "Kc_mid", "value": f"{kc_mid:.3f} (table {kc_mid_tab:.3f})",
                               "source": f"FAO-56 eq. 62 with mid-season mean u2 {u2m:.2f} m/s, RHmin {rhm:.0f}%, "
                                         f"h {h:g} m ({len(mid_days)} days)"})
            else:
                warnings.append("No RHmin/wind data in the mid-season: Kc_mid was not adjusted for climate.")
            u2l, rhl = stage_means(late_days)
            if kc_end_tab >= 0.45 and u2l is not None and rhl is not None:
                kc_end = kc_climate_adjust(kc_end_tab, u2l, rhl, h)
                params.append({"parameter": "Kc_end", "value": f"{kc_end:.3f} (table {kc_end_tab:.3f})",
                               "source": f"FAO-56 eq. 65 with late-season mean u2 {u2l:.2f} m/s, RHmin {rhl:.0f}%"})
        if not any(p["parameter"] == "Kc_mid" for p in params):
            params.append({"parameter": "Kc_mid", "value": f"{kc_mid:.3f}", "source": "user" if args.kc_mid is not None else src})
        if not any(p["parameter"] == "Kc_end" for p in params):
            why = " (not adjusted: below 0.45, eq. 65 applies only to crops harvested green)" if kc_end_tab < 0.45 else ""
            params.append({"parameter": "Kc_end", "value": f"{kc_end:.3f}",
                           "source": ("user" if args.kc_end is not None else src) + why})
        params += [{"parameter": "Crop height h (m)", "value": "" if h is None else f"{h:g}", "source": src},
                   {"parameter": "Stage lengths (days)", "value": ",".join(map(str, stages)), "source": stage_src},
                   {"parameter": "Root depth Zr (m)", "value": f"{zr:g}", "source": zr_src},
                   {"parameter": "p (ETc ~ 5 mm/day)", "value": f"{p_table:g}",
                    "source": "user" if args.p is not None else "FAO-56 Table 22"},
                   {"parameter": "TAW (mm)", "value": f"{taw:.1f}", "source": taw_src}]
        if crop.get("note"):
            params.append({"parameter": "Note", "value": "", "source": crop["note"]})

        if args.initial_depletion_mm == 0:
            warnings.append("Initial depletion was set to 0 (root zone at field capacity on the planting day). "
                            "Use --initial-depletion-mm if the spring was dry.")
        for w in days:
            i = (w["date"] - plant).days + 1
            w["kc"], w["stage"] = kc_on_day(i, stages, kc_ini, kc_mid, kc_end)
        if args.mad is not None and not 0.05 <= args.mad <= 0.95:
            raise ValueError("--mad must be a fraction of TAW between 0.05 and 0.95.")
        rows = run_balance(days, taw, p_table, args.initial_depletion_mm, p_adjust=not args.no_p_adjust,
                           auto_irrigate=args.irrigate, mad=args.mad,
                           min_rain_frac_et0=args.ignore_rain_below_et0_frac)
        for r, w in zip(rows, days):
            r["stage"] = w["stage"]
        resid = mass_balance_residual(rows, args.initial_depletion_mm)
        if abs(resid) > 1e-6:
            raise RuntimeError(f"Internal error: water balance does not close (residual {resid} mm).")
        if any(r["et_limited_by_storage"] for r in rows):
            warnings.append("On some days ETc adj was limited by the water left in the root zone (Dr reached TAW).")
        if any(r["irrigation_mm"] > 0 and r["auto_irrigation_mm"] == 0 for r in rows):
            warnings.append("Irrigation from the weather file's irrigation_mm column was added as net depth.")

        stress = [r["stressed"] for r in rows]
        first_stress = next((r["date"].isoformat() for r in rows if r["stressed"]), "")
        summary = {"days": len(rows), "first_date": rows[0]["date"].isoformat(), "last_date": rows[-1]["date"].isoformat(),
                   **{k: sum(r[k] for r in rows) for k in ("precip_mm", "runoff_mm", "precip_ignored_mm", "eff_precip_mm",
                                                           "irrigation_mm", "capillary_rise_mm", "etc_mm", "etc_adj_mm",
                                                           "deep_percolation_mm")},
                   "stress_days": sum(stress), "first_stress_date": first_stress, "longest_stress_run_days": runs(stress),
                   "irrigation_events": sum(1 for r in rows if r["auto_irrigation_mm"] > 0),
                   "final_depletion_mm": rows[-1]["dr_end_mm"], "taw_mm": taw, "mass_balance_residual_mm": resid}
        stress_by_stage = {n: sum(1 for r in rows if r["stage"] == n and r["stressed"]) for n in STAGE_NAMES}
        events = [{"date": r["date"].isoformat(), "net_irrigation_mm": r["auto_irrigation_mm"],
                   "depletion_before_mm": r["dr_start_mm"]} for r in rows if r["auto_irrigation_mm"] > 0]

        out = Path(args.out)
        out.mkdir(parents=True, exist_ok=True)
        table = [{"date": r["date"].isoformat(), "day": (r["date"] - plant).days + 1, "stage": r["stage"],
                  **{k: (round(v, 4) if isinstance(v, float) else ("yes" if v is True else "no" if v is False else v))
                     for k, v in r.items() if k not in ("date", "stage")}} for r in rows]
        save_table(table, out, "daily_water_balance.csv")
        save_table([{k: (round(v, 3) if isinstance(v, float) else v) for k, v in summary.items()}], out,
                   "season_summary.csv")
        save_table(params, out, "crop_parameters.csv")
        save_table(events or [{"date": "", "net_irrigation_mm": "", "depletion_before_mm": ""}], out,
                   "irrigation_events.csv")
        bounds = [("development", plant + timedelta(days=stages[0])),
                  ("mid-season", plant + timedelta(days=stages[0] + stages[1])),
                  ("late season", plant + timedelta(days=stages[0] + stages[1] + stages[2]))]
        figures = make_figures(rows, [b for b in bounds if b[1] <= rows[-1]["date"]], out)
        irr_desc = ("none" if not args.irrigate else
                    f"scheduled when depletion ≥ {'MAD ' + str(args.mad) if args.mad is not None else 'RAW (MAD = p)'} "
                    "× TAW; net depth refills to field capacity")
        write_report(out, {"warnings": warnings, "crop_name": crop.get("name", args.crop), "plant": plant.isoformat(),
                           "stages": stages, "season_end": season_end.isoformat(), "weather": args.weather,
                           "summary": summary, "zr": zr, "taw": taw, "taw_source": taw_src, "p_table": p_table,
                           "p_adjust": not args.no_p_adjust, "dr0": args.initial_depletion_mm,
                           "irrigation_desc": irr_desc, "params": params, "events": events,
                           "stress_by_stage": stress_by_stage, "figures": figures})
        write_repro_bundle(out, argv, inputs=inputs, packages=["matplotlib"])
    except (ValueError, RuntimeError) as e:
        print(ascii_safe(f"Error: {e}"), file=sys.stderr)
        return 1
    print(ascii_safe(f"Done: {summary['days']} days, ETc {summary['etc_mm']:.0f} mm, ETc adj "
                     f"{summary['etc_adj_mm']:.0f} mm, {summary['stress_days']} stress days. "
                     f"Report: {(out / 'report.md').resolve()}"))
    for w in warnings:
        print(ascii_safe(f"Warning: {w}"), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
