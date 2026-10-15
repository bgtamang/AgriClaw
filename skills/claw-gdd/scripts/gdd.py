#!/usr/bin/env python3
"""claw-gdd: growing degree days (GDD) from daily maximum and minimum air
temperature, corn stage timing from a cited table in its native units
(degF-days), daylength for photoperiod-sensitive crops, and freeze dates.

Usage:
    python scripts/gdd.py --demo
    python scripts/gdd.py --lat 40.08 --lon -88.24 --plant-date 2025-04-25 --crop corn
    python scripts/gdd.py --input station.csv --temp-units F --lat 40.08 \
        --plant-date 2025-04-25 --crop corn

Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT license.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import platform
import shlex
import sys
import urllib.error
import urllib.parse
import urllib.request
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

SKILL = "claw-gdd"
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
# Constants (sources in SKILL.md, section "Expert decisions")
# ---------------------------------------------------------------------------

# Crop presets. Corn: the US 86/50 degF cut-off ("modified") method of
# Neild & Newman (1987, NCH-40). Soybean: the same thresholds, used here only as
# a thermal-time convention (no stage prediction; see SKILL.md).
CROP_PRESETS = {
    "corn": {"base_c": 10.0, "upper_c": 30.0, "method": "modified"},
    "soybean": {"base_c": 10.0, "upper_c": 30.0, "method": "modified"},
}

# Neild & Newman (1987, Purdue reprint RR 4/90), NCH-40 Table 3: GDD from planting
# (degF-days, base 50 degF, 86/50 cut-off method) for a 2700 GDD hybrid.
# (Hanway stage number, leaf-collar / R code, NCH-40 description, degF-days)
CORN_STAGES_NCH40 = [
    (0.5, "V2", "Two leaves fully emerged", 200),
    (1.0, "V4", "Four leaves fully emerged", 345),
    (1.5, "V6", "Six leaves fully emerged (growing point above soil)", 475),
    (2.0, "V8", "Eight leaves fully emerged (tassel beginning to develop)", 610),
    (2.5, "V10", "Ten leaves fully emerged", 740),
    (3.0, "V12", "Twelve leaves fully emerged (ear formation)", 870),
    (3.5, "V14", "Fourteen leaves fully emerged (silks developing on ear)", 1000),
    (4.0, "V16", "Sixteen leaves fully emerged (tip of tassel emerging)", 1135),
    (5.0, "R1", "Silks emerging / pollen shedding (plant at full height)", 1400),
    (6.0, "R2", "Kernels in blister stage", 1660),
    (7.0, "R4", "Kernels in dough stage", 1925),
    (8.0, "R5", "Kernels denting", 2190),
    (9.0, "R5 (late)", "Kernels dented", 2450),
    (10.0, "R6", "Physiological maturity (black layer)", 2700),
]
CORN_TABLE_RATING_F = 2700.0

# NCEI frost-freeze thresholds are "at or below" 32 degF (freeze) and 28 degF
# (hard freeze); NCEI computes freeze dates over an August 1 - July 31 year.
FREEZE_C = 0.0
HARD_FREEZE_C = (28.0 - 32.0) / 1.8

MAX_GAP_DAYS = 3          # gaps of this many days or fewer are interpolated
DEFAULT_SEASON_DAYS = 210  # default end date = planting + this, clamped to data
PAST_DAYS_MAX = 92         # Open-Meteo forecast API limit for &past_days=

ARCHIVE_URL = "https://archive-api.open-meteo.com/v1/archive"
FORECAST_URL = "https://api.open-meteo.com/v1/forecast"
OM_VARS = {"temperature_2m_max": ("tmax_c", "°C"), "temperature_2m_min": ("tmin_c", "°C")}
MODELS = ("best_match", "era5", "era5_seamless", "era5_land", "ecmwf_ifs")

DEMO = {"lat": 40.08, "lon": -88.24, "plant_date": "2025-04-25", "crop": "corn",
        "file": "openmeteo_champaign_2025.json"}


class GapError(ValueError):
    """Raised when weather data have a gap that is too long to interpolate."""


# ---------------------------------------------------------------------------
# Units, GDD, daylength
# ---------------------------------------------------------------------------


def c_to_f(t_c: float) -> float:
    return t_c * 1.8 + 32.0


def f_to_c(t_f: float) -> float:
    return (t_f - 32.0) / 1.8


def daily_gdd(tmax: float, tmin: float, base: float, upper: float | None = None,
              method: str = "modified") -> float:
    """Daily growing degree days (McMaster & Wilhelm 1997).

    average:  max(0, (Tmax + Tmin)/2 - Tbase); no upper threshold.
    modified: Tmax and Tmin are each limited to [Tbase, Tupper] before averaging
              (the US corn 86/50 degF cut-off method). Without Tupper only the
              lower limit applies.
    Units follow the inputs: degC in -> degC-days, degF in -> degF-days.
    """
    if method == "average":
        return max(0.0, (tmax + tmin) / 2.0 - base)
    if method == "modified":
        hi = math.inf if upper is None else upper
        tx = min(max(tmax, base), hi)
        tn = min(max(tmin, base), hi)
        return (tx + tn) / 2.0 - base
    raise ValueError(f"Unknown GDD method '{method}'; use 'average' or 'modified'.")


def daylength_hours(lat_deg: float, doy: int, p_deg: float = 0.8333) -> float:
    """Daylength (h) from the CBM model of Forsythe et al. (1995).

    p_deg = 0.8333 deg: sunrise/sunset with the top of the sun at the horizon,
    including refraction; p_deg = 6: includes civil twilight.
    """
    theta = 0.2163108 + 2.0 * math.atan(0.9671396 * math.tan(0.00860 * (doy - 186)))
    phi = math.asin(0.39795 * math.cos(theta))
    lat = math.radians(lat_deg)
    den = math.cos(lat) * math.cos(phi)
    if abs(den) < 1e-12:
        return 24.0 if math.sin(lat) * math.sin(phi) > 0 else 0.0
    x = (math.sin(math.radians(p_deg)) + math.sin(lat) * math.sin(phi)) / den
    x = max(-1.0, min(1.0, x))
    return 24.0 - (24.0 / math.pi) * math.acos(x)


def parse_date(text: str, what: str = "date") -> date:
    try:
        return date.fromisoformat(str(text).strip())
    except ValueError:
        raise ValueError(f"{what} '{text}' is not a valid YYYY-MM-DD date.") from None


# ---------------------------------------------------------------------------
# Open-Meteo access (fetching is separate from parsing so tests run offline)
# ---------------------------------------------------------------------------


def openmeteo_url(kind: str, lat: float, lon: float, start: date | None = None,
                  end: date | None = None, model: str = "best_match",
                  past_days: int | None = None, variables=tuple(OM_VARS)) -> str:
    params = {"latitude": f"{lat:.5f}", "longitude": f"{lon:.5f}",
              "daily": ",".join(variables), "timezone": "auto"}
    if kind == "archive":
        params.update(start_date=start.isoformat(), end_date=end.isoformat())
        if model != "best_match":
            params["models"] = model
        return ARCHIVE_URL + "?" + urllib.parse.urlencode(params)
    if kind == "forecast":
        params.update(past_days=str(past_days), forecast_days="0")
        return FORECAST_URL + "?" + urllib.parse.urlencode(params)
    raise ValueError(kind)


def http_get_json(url: str, timeout: int = 90) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": f"AgriClaw/{VERSION} ({SKILL})"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", "replace")
        try:
            reason = json.loads(body).get("reason", body)
        except json.JSONDecodeError:
            reason = body[:300]
        raise RuntimeError(f"Open-Meteo returned HTTP {e.code}: {reason}") from None
    except (urllib.error.URLError, TimeoutError) as e:
        raise RuntimeError(f"Could not reach Open-Meteo ({e}). Check the internet "
                           "connection, or use --input with your own station data.") from None


def parse_openmeteo_daily(resp: dict, source: str) -> tuple[list[dict], dict]:
    """Rows of {date, tmax_c, tmin_c, source} plus response metadata."""
    if resp.get("error"):
        raise RuntimeError(f"Open-Meteo error: {resp.get('reason', resp)}")
    daily = resp.get("daily")
    if not isinstance(daily, dict) or "time" not in daily:
        raise RuntimeError("Open-Meteo response has no 'daily' block; cannot parse it.")
    units = resp.get("daily_units", {})
    for var, (_, unit) in OM_VARS.items():
        if var not in daily:
            raise RuntimeError(f"Open-Meteo response lacks '{var}'.")
        if units.get(var) not in (None, unit):
            raise RuntimeError(f"'{var}' is in {units.get(var)}, expected {unit}.")
    rows = []
    for i, t in enumerate(daily["time"]):
        row = {"date": parse_date(t)}
        for var, (col, _) in OM_VARS.items():
            v = daily[var][i]
            row[col] = None if v is None else float(v)
        row["source"] = source
        rows.append(row)
    meta = {k: resp.get(k) for k in ("latitude", "longitude", "elevation", "timezone",
                                     "utc_offset_seconds")}
    return rows, meta


def local_today(utc_offset_seconds: int | None, now_utc: datetime | None = None) -> date:
    now_utc = now_utc or datetime.now(timezone.utc)
    return (now_utc + timedelta(seconds=utc_offset_seconds or 0)).date()


def trim_trailing_missing(rows: list[dict], cols=("tmax_c", "tmin_c")) -> tuple[list[dict], int]:
    """Drop trailing days where every column is missing (data not yet published)."""
    n = len(rows)
    while n and all(rows[n - 1][c] is None for c in cols):
        n -= 1
    return rows[:n], len(rows) - n


def fetch_daily(lat: float, lon: float, start: date, end: date, model: str = "best_match",
                today: date | None = None, fetcher=http_get_json,
                now_utc: datetime | None = None) -> tuple[list[dict], dict, list[str], list[dict]]:
    """Fetch daily Tmax/Tmin for [start, end] without ever asking for future days.

    The end date is clamped to yesterday. The archive API is queried first; if its
    most recent days are not yet available (ERA5 is published with a delay of about
    5 days), the gap up to yesterday is filled from the forecast API with
    &past_days= (at most 92 days back). Returns rows, metadata, notes and the raw
    responses (for saving in the output folder).
    """
    today = today or date.today()
    notes: list[str] = []
    raw: list[dict] = []
    yesterday = today - timedelta(days=1)
    if start > yesterday:
        raise ValueError(f"Start date {start} is not in the past (today is {today}); "
                         "there are no observations to use yet.")
    end_eff = min(end, yesterday)
    if end_eff < end:
        notes.append(f"Requested end date {end} is after the last complete day; the "
                     f"analysis stops at {end_eff} at the latest. Nothing is forecast "
                     "or extrapolated.")
    url = openmeteo_url("archive", lat, lon, start, end_eff, model=model)
    resp = fetcher(url)
    raw.append({"url": url, "response": resp})
    rows, meta = parse_openmeteo_daily(resp, f"archive:{model}")
    ltoday = local_today(meta.get("utc_offset_seconds"), now_utc)
    rows = [r for r in rows if r["date"] < ltoday]
    rows, n_trim = trim_trailing_missing(rows)
    if n_trim:
        notes.append(f"The archive ({model}) has no data yet for the last {n_trim} day(s) "
                     "requested (reanalysis is published with a delay).")
    last = rows[-1]["date"] if rows else start - timedelta(days=1)
    if last < end_eff:
        first_missing = last + timedelta(days=1)
        past_days = (today - first_missing).days
        if past_days <= PAST_DAYS_MAX:
            url2 = openmeteo_url("forecast", lat, lon, past_days=past_days)
            resp2 = fetcher(url2)
            raw.append({"url": url2, "response": resp2})
            rows2, _ = parse_openmeteo_daily(resp2, "forecast:past_days")
            rows2 = [r for r in rows2 if first_missing <= r["date"] <= end_eff and r["date"] < ltoday]
            rows2, _ = trim_trailing_missing(rows2)
            if rows2:
                rows += rows2
                notes.append(f"{len(rows2)} most recent day(s) ({rows2[0]['date']} to "
                             f"{rows2[-1]['date']}) come from the Open-Meteo forecast API "
                             "(&past_days=): operational model output, provisional, and not "
                             "the same dataset as the archive.")
        else:
            notes.append(f"Data after {last} are not available from the archive yet and are "
                         "too old for the forecast API's past_days window.")
    return rows, meta, notes, raw


# ---------------------------------------------------------------------------
# Station CSV input
# ---------------------------------------------------------------------------

MISSING_TOKENS = {"", "na", "nan", "null", "none", "missing"}


def read_station_csv(path: str, temp_units: str) -> list[dict]:
    """Read date, tmax, tmin from a CSV. Temperature units must be given
    explicitly (temp_units 'C' or 'F'); column suffixes _c/_f must agree."""
    p = Path(path)
    if not p.is_file():
        raise ValueError(f"Input file not found: {path}")
    if temp_units not in ("C", "F"):
        raise ValueError("--temp-units must be C or F.")
    with p.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        header = {h.strip().lower(): h for h in (reader.fieldnames or [])}
        raw_rows = list(reader)
    if "date" not in header:
        raise ValueError(f"No 'date' column in {path}. Columns found: {', '.join(header) or 'none'}.")
    cols = {}
    for key in ("tmax", "tmin"):
        found = [h for h in (key, f"{key}_c", f"{key}_f") if h in header]
        if len(found) != 1:
            raise ValueError(f"Need exactly one of {key}, {key}_c or {key}_f in {path}; found "
                             f"{found or 'none'}.")
        h = found[0]
        if h.endswith("_c") and temp_units == "F" or h.endswith("_f") and temp_units == "C":
            raise ValueError(f"Column '{header[h]}' says {h[-1].upper()} but --temp-units is "
                             f"{temp_units}. Check which unit the data are in.")
        cols[key] = header[h]
    rows, seen = [], set()
    for i, r in enumerate(raw_rows, start=2):
        d = parse_date(r[header["date"]], f"Row {i}: date")
        if d in seen:
            raise ValueError(f"Row {i}: date {d} appears more than once.")
        seen.add(d)
        row = {"date": d, "source": "input CSV"}
        for key in ("tmax", "tmin"):
            raw_v = (r[cols[key]] or "").strip()
            if raw_v.lower() in MISSING_TOKENS:
                row[f"{key}_c"] = None
                continue
            try:
                v = float(raw_v)
            except ValueError:
                raise ValueError(f"Row {i}: {cols[key]} '{raw_v}' is not numeric.") from None
            row[f"{key}_c"] = f_to_c(v) if temp_units == "F" else v
        rows.append(row)
    rows.sort(key=lambda r: r["date"])
    check_temperature_sanity(rows)
    return rows


def check_temperature_sanity(rows: list[dict]) -> None:
    bad = [r for r in rows if any(r[c] is not None and not -60.0 <= r[c] <= 60.0
                                  for c in ("tmax_c", "tmin_c"))]
    if bad:
        r = bad[0]
        raise ValueError(f"{len(bad)} day(s) have temperatures outside -60..60 degC (first: {r['date']}, "
                         f"tmax {r['tmax_c']}, tmin {r['tmin_c']}). Check --temp-units, and "
                         "replace missing-value codes such as -9999 with blank cells.")
    swapped = [r for r in rows if r["tmax_c"] is not None and r["tmin_c"] is not None
               and r["tmin_c"] > r["tmax_c"] + 0.05]
    if swapped:
        raise ValueError(f"{len(swapped)} day(s) have tmin > tmax (first: {swapped[0]['date']}). "
                         "Check that the tmax and tmin columns are not swapped.")


# ---------------------------------------------------------------------------
# Calendar completion and gap handling
# ---------------------------------------------------------------------------


def complete_calendar(rows: list[dict], start: date, end: date) -> list[dict]:
    """One row per day from start to end; days absent from the data get None."""
    by_date = {r["date"]: r for r in rows}
    out, d = [], start
    while d <= end:
        out.append(dict(by_date.get(d) or {"date": d, "tmax_c": None, "tmin_c": None,
                                            "source": "missing"}))
        d += timedelta(days=1)
    return out


def fill_short_gaps(rows: list[dict], cols=("tmax_c", "tmin_c"),
                    max_gap: int = MAX_GAP_DAYS) -> tuple[list[dict], list[dict], list[dict]]:
    """Linearly interpolate runs of missing values of <= max_gap days that have data
    on both sides. Returns (rows, filled, unfilled_gaps); nothing is filled silently."""
    rows = [dict(r, filled="") for r in rows]
    filled, gaps = [], []
    for c in cols:
        i, n = 0, len(rows)
        while i < n:
            if rows[i][c] is not None:
                i += 1
                continue
            j = i
            while j < n and rows[j][c] is None:
                j += 1
            length = j - i
            if 0 < i and j < n and length <= max_gap:
                a, b = rows[i - 1][c], rows[j][c]
                for k in range(i, j):
                    frac = (k - i + 1) / (length + 1)
                    rows[k][c] = a + frac * (b - a)
                    rows[k]["filled"] = (rows[k]["filled"] + " " + c).strip()
                    filled.append({"date": rows[k]["date"], "variable": c,
                                   "value_c": round(rows[k][c], 2), "gap_days": length})
            else:
                gaps.append({"variable": c, "start": rows[i]["date"], "end": rows[j - 1]["date"],
                             "days": length,
                             "reason": "longer than max gap" if length > max_gap else "at edge of data"})
            i = j
    return rows, filled, gaps


# ---------------------------------------------------------------------------
# Season analysis
# ---------------------------------------------------------------------------


def accumulate_gdd(rows: list[dict], plant: date, end: date, base_c: float,
                   upper_c: float | None, method: str, lat: float) -> list[dict]:
    """Daily table from planting (day 1, inclusive) to end."""
    out, cum = [], 0.0
    for r in rows:
        if not plant <= r["date"] <= end:
            continue
        if r["tmax_c"] is None or r["tmin_c"] is None:
            raise GapError(f"Missing temperature on {r['date']} inside the season.")
        g = daily_gdd(r["tmax_c"], r["tmin_c"], base_c, upper_c, method)
        cum += g
        doy = r["date"].timetuple().tm_yday
        out.append({
            "date": r["date"].isoformat(),
            "day_after_planting": (r["date"] - plant).days + 1,
            "tmax_c": round(r["tmax_c"], 2), "tmin_c": round(r["tmin_c"], 2),
            "tmax_f": round(c_to_f(r["tmax_c"]), 2), "tmin_f": round(c_to_f(r["tmin_c"]), 2),
            "gdd_c": round(g, 3), "gdd_f": round(g * 1.8, 3),
            "cum_gdd_c": round(cum, 2), "cum_gdd_f": round(cum * 1.8, 2),
            "daylength_h": round(daylength_hours(lat, doy), 3),
            "daylength_civil_twilight_h": round(daylength_hours(lat, doy, 6.0), 3),
            "filled": r.get("filled", ""), "source": r.get("source", ""),
        })
    return out


def corn_stages(daily: list[dict], hybrid_rating_f: float | None = None) -> list[dict]:
    """First date on which cumulative degF-days reach each NCH-40 threshold.
    Thresholds are scaled by hybrid_rating_f / 2700 when a rating is given."""
    scale = (hybrid_rating_f / CORN_TABLE_RATING_F) if hybrid_rating_f else 1.0
    last = daily[-1] if daily else None
    out = []
    for num, code, desc, thr_f in CORN_STAGES_NCH40:
        thr = thr_f * scale
        hit = next((d for d in daily if d["cum_gdd_f"] >= thr - 1e-9), None)
        out.append({
            "hanway_stage": num, "stage": code, "description_nch40": desc,
            "gdd_f_table": thr_f, "gdd_f_used": round(thr, 1), "gdd_c_used": round(thr / 1.8, 1),
            "date_reached": hit["date"] if hit else "",
            "days_after_planting": hit["day_after_planting"] if hit else "",
            "status": "reached" if hit else
            (f"not reached as of {last['date']} ({last['cum_gdd_f']:.0f} of {thr:.0f} degF-days)"
             if last else "no data"),
        })
    return out


def fall_boundary(plant: date, lat: float) -> date:
    """Start of the 'fall' frost window: August 1 (northern hemisphere) or
    February 1 (southern) after planting, mirroring NCEI's Aug 1 - Jul 31 frost year."""
    if lat >= 0:
        return date(plant.year if plant < date(plant.year, 8, 1) else plant.year + 1, 8, 1)
    return date(plant.year + (1 if plant >= date(plant.year, 2, 1) else 0), 2, 1)


def frost_analysis(rows: list[dict], plant: date, lat: float, freeze_c: float = FREEZE_C,
                   hard_c: float = HARD_FREEZE_C) -> tuple[list[dict], dict]:
    """Freeze events (Tmin <= threshold) and season summary dates."""
    boundary = fall_boundary(plant, lat)
    events = []
    for r in rows:
        t = r.get("tmin_c")
        if t is None or t > freeze_c:
            continue
        d = r["date"]
        period = ("before_planting" if d < plant else
                  "after_planting_before_fall" if d < boundary else "fall")
        events.append({"date": d.isoformat(), "tmin_c": round(t, 2), "tmin_f": round(c_to_f(t), 1),
                       "hard_freeze": "yes" if t <= hard_c else "no", "period": period})

    def pick(period, hard, last):
        ev = [e for e in events if e["period"] == period and (not hard or e["hard_freeze"] == "yes")]
        if not ev:
            return ""
        return ev[-1]["date"] if last else ev[0]["date"]

    valid = [r["date"] for r in rows if r.get("tmin_c") is not None]
    summary = {
        "search_start": valid[0].isoformat() if valid else "",
        "search_end": valid[-1].isoformat() if valid else "",
        "fall_window_starts": boundary.isoformat(),
        "last_freeze_before_planting": pick("before_planting", False, True),
        "last_hard_freeze_before_planting": pick("before_planting", True, True),
        "freezes_after_planting_before_fall": sum(e["period"] == "after_planting_before_fall"
                                                  for e in events),
        "last_freeze_after_planting_before_fall": pick("after_planting_before_fall", False, True),
        "first_fall_freeze": pick("fall", False, False),
        "first_fall_hard_freeze": pick("fall", True, False),
        "fall_window_observed": bool(valid) and valid[-1] >= boundary,
    }
    return events, summary


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def format_latlon(lat: float | None, lon: float | None) -> str:
    if lat is None or lon is None:
        return "not given"
    return (f"{abs(lat):.4f}°{'N' if lat >= 0 else 'S'}, "
            f"{abs(lon):.4f}°{'E' if lon >= 0 else 'W'}")


def make_figures(daily: list[dict], stages: list[dict] | None, base_c: float,
                 upper_c: float | None, freeze_c: float, out: Path) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    fig_dir = out / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    dates = [date.fromisoformat(d["date"]) for d in daily]
    files = []

    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.plot(dates, [d["cum_gdd_f"] for d in daily], color="#2f6b3a", lw=2)
    ax.set_ylabel("Cumulative GDD (°F-days)")
    sec = ax.secondary_yaxis("right", functions=(lambda x: x / 1.8, lambda x: x * 1.8))
    sec.set_ylabel("Cumulative GDD (°C-days)")
    if stages:
        for s in stages:
            ax.axhline(s["gdd_f_used"], color="#999999", lw=0.6, ls=":")
            if s["date_reached"]:
                dd = date.fromisoformat(s["date_reached"])
                ax.plot([dd], [s["gdd_f_used"]], "o", color="#b35900", ms=4)
                ax.annotate(s["stage"], (dd, s["gdd_f_used"]), textcoords="offset points",
                            xytext=(-4, 4), ha="right", fontsize=7)
    ax.set_title("Growing degree-day accumulation from planting")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax.grid(alpha=0.3)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(fig_dir / "gdd_accumulation.png", dpi=150)
    plt.close(fig)
    files.append("figures/gdd_accumulation.png")

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9, 6), sharex=True,
                                 gridspec_kw={"height_ratios": [3, 1.3]})
    a1.plot(dates, [d["tmax_c"] for d in daily], color="#c0392b", lw=1, label="Tmax")
    a1.plot(dates, [d["tmin_c"] for d in daily], color="#2e6da4", lw=1, label="Tmin")
    a1.axhline(base_c, color="#555555", ls="--", lw=0.8, label=f"Tbase {base_c:g} °C")
    if upper_c is not None:
        a1.axhline(upper_c, color="#555555", ls="-.", lw=0.8, label=f"Tupper {upper_c:g} °C")
    a1.axhline(freeze_c, color="#7f3fbf", ls=":", lw=0.8, label=f"freeze {freeze_c:g} °C")
    filled = [(dd, d["tmax_c"]) for dd, d in zip(dates, daily) if d["filled"]]
    if filled:
        a1.plot([f[0] for f in filled], [f[1] for f in filled], "x", color="black",
                ms=5, label="interpolated")
    a1.set_ylabel("Air temperature (°C)")
    a1.legend(fontsize=7, ncol=3, loc="lower center")
    a1.grid(alpha=0.3)
    a2.plot(dates, [d["daylength_h"] for d in daily], color="#8a6d00", lw=1.2, label="sunrise-sunset")
    a2.plot(dates, [d["daylength_civil_twilight_h"] for d in daily], color="#8a6d00", lw=1,
            ls="--", label="incl. civil twilight")
    a2.set_ylabel("Daylength (h)")
    a2.legend(fontsize=7)
    a2.grid(alpha=0.3)
    a2.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(fig_dir / "temperature_daylength.png", dpi=150)
    plt.close(fig)
    files.append("figures/temperature_daylength.png")
    return files


REFERENCES = [
    "McMaster, G.S. & Wilhelm, W.W. (1997). Growing degree-days: one equation, two interpretations. "
    "Agricultural and Forest Meteorology 87(4), 291-300. https://doi.org/10.1016/S0168-1923(97)00027-0",
    "Neild, R.E. & Newman, J.E. (1987; Purdue reprint RR 4/90). Growing Season Characteristics and "
    "Requirements in the Corn Belt. National Corn Handbook NCH-40, Purdue University Cooperative "
    "Extension Service. https://www.extension.purdue.edu/extmedia/nch/nch-40.html",
    "Lee, C. (2011). Corn Growth Stages and Growing Degree Days: A Quick Reference Guide. AGR-202, "
    "University of Kentucky Cooperative Extension Service. "
    "https://publications.ca.uky.edu/sites/publications.ca.uky.edu/files/agr202.pdf",
    "Forsythe, W.C., Rykiel, E.J., Stahl, R.S., Wu, H. & Schoolfield, R.M. (1995). A model comparison "
    "for daylength as a function of latitude and day of year. Ecological Modelling 80(1), 87-95. "
    "https://doi.org/10.1016/0304-3800(94)00034-F",
    "Arguez, A. (2012). NOAA's 1981-2010 Climate Normals: Supplemental Normals. NOAA National Climatic "
    "Data Center technical report (frost-freeze thresholds and the August 1 - July 31 frost year), "
    "reproduced in Appendix I of NOAA NCEI, Normals Calculation Methodology 2020. https://www.ncei.noaa.gov/"
    "data/normals-annualseasonal/1991-2020/doc/Normals_Calculation_Methodology_2020.pdf",
    "Setiyono, T.D., Weiss, A., Specht, J., Bastidas, A.M., Cassman, K.G. & Dobermann, A. (2007). "
    "Understanding and modeling the effect of temperature and daylength on soybean phenology under "
    "high-yield conditions. Field Crops Research 100(2-3), 257-271. https://doi.org/10.1016/j.fcr.2006.07.011",
    "Salmerón, M. & Purcell, L.C. (2016). Simplifying the prediction of phenology with the "
    "DSSAT-CROPGRO-soybean model based on relative maturity group and determinacy. Agricultural "
    "Systems 148, 178-187. https://doi.org/10.1016/j.agsy.2016.07.016",
    "Zippenfenig, P. (2023). Open-Meteo.com Weather API [Computer software]. Zenodo. "
    "https://doi.org/10.5281/zenodo.7970649",
    "Hersbach, H. et al. (2020). The ERA5 global reanalysis. Quarterly Journal of the Royal "
    "Meteorological Society 146(730), 1999-2049. https://doi.org/10.1002/qj.3803",
]


def write_report(out: Path, ctx: dict) -> None:
    L = [f"# {SKILL} report", ""]
    if ctx["warnings"]:
        L += ["## Warnings", ""] + [f"- {w}" for w in ctx["warnings"]] + [""]
    L += ["## Inputs", "",
          f"- Location: {format_latlon(ctx['lat'], ctx['lon'])}",
          f"- Weather source: {ctx['source_desc']}",
          f"- Crop: {ctx['crop']}; planting date {ctx['plant']} (counted as day 1)",
          f"- Season analyzed: {ctx['plant']} to {ctx['end']} ({ctx['n_days']} days)",
          f"- Method: {ctx['method']}; Tbase {ctx['base_c']:g} °C ({c_to_f(ctx['base_c']):g} °F)"
          + (f"; Tupper {ctx['upper_c']:g} °C ({c_to_f(ctx['upper_c']):g} °F)"
             if ctx["upper_c"] is not None and ctx["method"] == "modified" else "; no upper threshold"),
          "", "## Season summary", "",
          "| Metric | Value |", "|---|---|",
          f"| Cumulative GDD | {ctx['cum_c']:.0f} °C-days = {ctx['cum_f']:.0f} °F-days |",
          f"| Mean Tmax / Tmin | {ctx['mean_tmax']:.1f} / {ctx['mean_tmin']:.1f} °C |",
          f"| Days interpolated | {ctx['n_filled_days']} |", ""]
    if ctx["stages"]:
        scale_txt = (f" scaled to a {ctx['hybrid']:.0f} °F-day hybrid (factor "
                     f"{ctx['hybrid'] / CORN_TABLE_RATING_F:.3f})" if ctx["hybrid"] else
                     " for a 2700 °F-day hybrid (unscaled)")
        L += ["## Corn stage timing", "",
              f"Thresholds from Neild & Newman (1987) NCH-40 Table 3 (°F-days, 86/50 °F method){scale_txt}. "
              "Dates are the first day the cumulative °F-days reach the threshold; nothing is "
              "extrapolated beyond the last day of data.", "",
              "| Stage | NCH-40 description | °F-days | Date | Days after planting |",
              "|---|---|---|---|---|"]
        for s in ctx["stages"]:
            L.append(f"| {s['stage']} | {s['description_nch40']} | {s['gdd_f_used']:.0f} | "
                     f"{s['date_reached'] or s['status']} | {s['days_after_planting']} |")
        L.append("")
    if ctx["crop"] == "soybean":
        L += ["## Soybean development", "",
              "No soybean V- or R-stage dates are predicted. Soybean flowering and maturity depend on "
              "photoperiod and maturity group as well as temperature (Setiyono et al. 2007; Salmerón & "
              "Purcell 2016), so a GDD threshold alone would give wrong dates. Use a photothermal model "
              "instead, such as SoyStage (University of Arkansas; https://soystage.uark.edu/, maturity "
              "groups 3-6, US Midsouth) or CROPGRO-Soybean in DSSAT. Daily daylength is in "
              "`tables/daily_gdd.csv` for use with such models.", ""]
    fs = ctx["frost"]
    L += ["## Freeze dates", "",
          f"Freeze = Tmin ≤ {ctx['freeze_c']:.2f} °C; hard freeze = Tmin ≤ {ctx['hard_c']:.2f} °C "
          "(NCEI 32 °F and 28 °F thresholds). Searched "
          f"{fs['search_start']} to {fs['search_end']}; the fall window starts {fs['fall_window_starts']}.", "",
          "| Item | Date |", "|---|---|",
          f"| Last freeze before planting | {fs['last_freeze_before_planting'] or 'none in searched period'} |",
          f"| Last hard freeze before planting | {fs['last_hard_freeze_before_planting'] or 'none in searched period'} |",
          f"| Freezes after planting, before the fall window | {fs['freezes_after_planting_before_fall']} "
          f"(last: {fs['last_freeze_after_planting_before_fall'] or 'none'}) |",
          f"| First fall freeze | {fs['first_fall_freeze'] or ('none observed through ' + fs['search_end'] if fs['fall_window_observed'] else 'fall window not reached in the data')} |",
          f"| First fall hard freeze | {fs['first_fall_hard_freeze'] or ('none observed through ' + fs['search_end'] if fs['fall_window_observed'] else 'fall window not reached in the data')} |",
          ""]
    L += ["## Figures", ""] + [f"![{f}]({f})" for f in ctx["figures"]] + [""]
    L += ["## Method notes and limitations", "",
          "- GDD use daily Tmax and Tmin only; they ignore water stress, soil temperature at seed depth "
          "and hourly temperature patterns.",
          "- °F-days = 1.8 × °C-days exactly when the thresholds are equivalent (50/86 °F = 10/30 °C).",
          "- Gridded weather (Open-Meteo) is model/reanalysis output for a grid cell, not a field station; "
          "frost on a clear night in a low spot can be colder than the grid value.",
          "- Gaps of up to 3 days are filled by linear interpolation and listed in `tables/filled_days.csv`; "
          "longer gaps inside the season stop the run.", ""]
    L += ["## References", ""] + [f"- {r}" for r in REFERENCES] + ["", f"*{SKILL} {VERSION}*"]
    write_text(out / "report.md", "\n".join(L) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="gdd.py", description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--demo", action="store_true",
                    help="corn planted 2025-04-25 near Champaign, IL (bundled Open-Meteo data)")
    ap.add_argument("--lat", type=float, help="latitude, decimal degrees (needed for daylength and frost season)")
    ap.add_argument("--lon", type=float, help="longitude, decimal degrees (needed to fetch weather)")
    ap.add_argument("--input", help="station CSV with date, tmax, tmin (use instead of fetching)")
    ap.add_argument("--temp-units", choices=["C", "F"], help="units of tmax/tmin in --input (required with --input)")
    ap.add_argument("--plant-date", help="planting date YYYY-MM-DD (required)")
    ap.add_argument("--end-date", help=f"last day to analyze (default: planting + {DEFAULT_SEASON_DAYS} days, "
                                       "clamped to the last day with data)")
    ap.add_argument("--crop", choices=["corn", "soybean", "other"], help="crop (required)")
    ap.add_argument("--method", choices=["average", "modified"], help="GDD method (default: modified for corn/soybean)")
    ap.add_argument("--base-c", type=float, help="base temperature, °C (required for --crop other)")
    ap.add_argument("--upper-c", type=float, help="upper threshold, °C (modified method only)")
    ap.add_argument("--hybrid-gdd", type=float,
                    help="corn hybrid GDD to black layer, °F-days (86/50 method); scales stage thresholds")
    ap.add_argument("--model", choices=MODELS, default="best_match", help="Open-Meteo archive model")
    ap.add_argument("--freeze-c", type=float, default=FREEZE_C, help="freeze threshold, °C (default 0 = 32 °F)")
    ap.add_argument("--hard-freeze-c", type=float, default=HARD_FREEZE_C,
                    help="hard-freeze threshold, °C (default -2.22 = 28 °F)")
    ap.add_argument("--out", default=f"agriclaw_output/{SKILL}", help="output folder")
    return ap


def main(argv: list[str] | None = None, today: date | None = None, fetcher=http_get_json) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = build_parser()
    args = ap.parse_args(argv)
    warnings: list[str] = []
    try:
        if args.demo:
            args.lat = args.lat if args.lat is not None else DEMO["lat"]
            args.lon = args.lon if args.lon is not None else DEMO["lon"]
            args.plant_date = args.plant_date or DEMO["plant_date"]
            args.crop = args.crop or DEMO["crop"]
        if not args.plant_date:
            ap.error("--plant-date is required (never guess a planting date; ask the user)")
        if not args.crop:
            ap.error("--crop is required: corn, soybean or other")
        if args.lat is None:
            ap.error("--lat is required (daylength and the frost season depend on it)")
        if not -90 <= args.lat <= 90 or (args.lon is not None and not -180 <= args.lon <= 180):
            ap.error("--lat must be within -90..90 and --lon within -180..180")
        if args.input and not args.temp_units:
            ap.error("--temp-units C or F is required with --input (units are never guessed)")
        if not (args.demo or args.input) and args.lon is None:
            ap.error("give --lat and --lon to fetch weather, or --input FILE.csv, or --demo")

        plant = parse_date(args.plant_date, "--plant-date")
        end = parse_date(args.end_date, "--end-date") if args.end_date else \
            plant + timedelta(days=DEFAULT_SEASON_DAYS)
        if end < plant:
            raise ValueError("--end-date is before --plant-date.")

        preset = CROP_PRESETS.get(args.crop, {})
        base_c = args.base_c if args.base_c is not None else preset.get("base_c")
        if base_c is None:
            raise ValueError("--base-c is required for --crop other: ask which base temperature "
                             "the user's phenology model or reference uses.")
        method = args.method or preset.get("method") or ("modified" if args.upper_c is not None else "average")
        upper_c = args.upper_c if args.upper_c is not None else preset.get("upper_c")
        if method == "average" and args.upper_c is not None:
            warnings.append("--upper-c is ignored by the 'average' method.")
        if method == "average":
            upper_c = None
        if upper_c is not None and upper_c <= base_c:
            raise ValueError("--upper-c must be above --base-c.")
        overridden = (args.base_c is not None or args.upper_c is not None or args.method is not None)

        frost_start = date(plant.year, 1, 1)
        raw_responses: list[dict] = []
        out = Path(args.out)
        out.mkdir(parents=True, exist_ok=True)
        inputs = []
        if args.input:
            rows = read_station_csv(args.input, args.temp_units)
            inputs.append(args.input)
            source_desc = f"station CSV `{args.input}` (temperatures in °{args.temp_units})"
            if not rows or rows[0]["date"] > plant:
                raise ValueError(f"The CSV starts on {rows[0]['date'] if rows else 'no date'}, after the "
                                 f"planting date {plant}; GDD from planting cannot be computed.")
            if rows[0]["date"] > frost_start:
                warnings.append(f"The CSV starts on {rows[0]['date']}, so the last-spring-freeze search "
                                f"covers {rows[0]['date']} to planting, not January 1.")
            data_end = rows[-1]["date"]
            rows = [r for r in rows if frost_start <= r["date"] <= end]
        elif args.demo:
            demo_file = SKILL_DIR / "demo" / DEMO["file"]
            resp = json.loads(demo_file.read_text(encoding="utf-8"))
            rows, meta = parse_openmeteo_daily(resp, "archive:best_match (bundled demo file)")
            inputs.append(str(demo_file))
            source_desc = (f"Open-Meteo archive (best match: ECMWF IFS 9 km / ERA5 / ERA5-Land), grid cell "
                           f"{format_latlon(meta['latitude'], meta['longitude'])}, saved in demo/{DEMO['file']}")
            rows = [r for r in rows if frost_start <= r["date"] <= end]
            data_end = rows[-1]["date"]
        else:
            if plant >= (today or date.today()):
                raise ValueError(f"Planting date {plant} is not in the past; there is no weather to "
                                 "accumulate yet (this skill never uses forecasts or climatology).")
            rows, meta, notes, raw_responses = fetch_daily(args.lat, args.lon, frost_start, end,
                                                           model=args.model, today=today, fetcher=fetcher)
            warnings += notes
            source_desc = (f"Open-Meteo archive, model '{args.model}' (gridded model/reanalysis data, not a "
                           f"station), grid cell {format_latlon(meta['latitude'], meta['longitude'])}")
            write_text(out / "raw" / "openmeteo_responses.json", json.dumps(raw_responses, indent=1))
            data_end = rows[-1]["date"] if rows else frost_start
        if not rows:
            raise ValueError("No weather data in the requested window.")

        rows = complete_calendar(rows, frost_start if not args.input else max(frost_start, rows[0]["date"]),
                                 min(end, data_end))
        rows, filled, gaps = fill_short_gaps(rows)
        # the season ends on the last day with both temperatures
        complete = [r["date"] for r in rows if r["tmax_c"] is not None and r["tmin_c"] is not None]
        end_eff = min(end, complete[-1]) if complete else plant
        if end_eff < end and not any("after the last complete day" in w for w in warnings):
            warnings.append(f"Data end on {end_eff}; the season was requested to {end}. Stages not reached "
                            f"by {end_eff} are reported as not reached (no extrapolation).")
        if end_eff < plant:
            raise ValueError(f"No complete weather data on or after the planting date {plant}.")
        season_gaps = [g for g in gaps if g["end"] >= plant and g["start"] <= end_eff]
        if season_gaps:
            g = season_gaps[0]
            raise GapError(f"{len(season_gaps)} gap(s) in {', '.join(sorted({x['variable'] for x in season_gaps}))} "
                           f"inside the season cannot be filled (first: {g['start']} to {g['end']}, {g['days']} "
                           f"days, {g['reason']}). This skill interpolates gaps of {MAX_GAP_DAYS} days or "
                           "fewer only: fill the gap from a nearby station or use the Open-Meteo data.")
        pre_gaps = [g for g in gaps if g["end"] < plant]
        if pre_gaps:
            warnings.append(f"{len(pre_gaps)} unfilled gap(s) before planting: the freeze search before "
                            "planting is incomplete.")
        if filled:
            warnings.append(f"{len({f['date'] for f in filled})} day(s) were filled by linear interpolation "
                            f"(gaps of {MAX_GAP_DAYS} days or fewer); see tables/filled_days.csv.")

        daily = accumulate_gdd(rows, plant, end_eff, base_c, upper_c, method, args.lat)
        stages = None
        if args.crop == "corn":
            standard = (method == "modified" and abs(base_c - 10.0) < 1e-9 and upper_c is not None
                        and abs(upper_c - 30.0) < 1e-9)
            if standard:
                if args.hybrid_gdd is not None:
                    if not 1800 <= args.hybrid_gdd <= 3600:
                        raise ValueError(f"--hybrid-gdd {args.hybrid_gdd:g} is outside 1800-3600 °F-days. "
                                         "Seed-company ratings are °F-days (86/50 method); a °C-day rating "
                                         "must be multiplied by 1.8, and Canadian Crop Heat Units (CHU) are "
                                         "a different scale that cannot be used here.")
                    warnings.append(f"Stage thresholds were scaled by {args.hybrid_gdd:g}/2700: proportional "
                                    "scaling of the NCH-40 table is a rule of thumb, not a calibrated model.")
                stages = corn_stages(daily, args.hybrid_gdd)
            else:
                warnings.append("Corn stage dates were not estimated: the NCH-40 thresholds are defined for the "
                                "86/50 °F modified method (Tbase 10 °C, Tupper 30 °C) and this run used "
                                f"{method}, Tbase {base_c:g} °C, Tupper {upper_c}.")
        elif args.hybrid_gdd is not None:
            warnings.append("--hybrid-gdd applies to corn only and was ignored.")
        if args.crop == "soybean":
            warnings.append("Soybean R stages are not predicted from GDD (photoperiod and maturity group "
                            "control them); see the 'Soybean development' section.")
        if overridden and args.crop in CROP_PRESETS:
            warnings.append(f"Default thresholds for {args.crop} were overridden by the user "
                            f"(method {method}, Tbase {base_c:g} °C, Tupper {upper_c}).")

        events, frost = frost_analysis([r for r in rows if r["date"] <= end_eff], plant, args.lat,
                                       args.freeze_c, args.hard_freeze_c)
        if args.lat < 0:
            warnings.append("Southern-hemisphere site: the fall freeze window starts on February 1 after "
                            "planting (mirror of NCEI's August 1 boundary).")

        save_table(daily, out, "daily_gdd.csv")
        if stages:
            save_table(stages, out, "corn_stages.csv")
        save_table(events or [{"date": "", "tmin_c": "", "tmin_f": "", "hard_freeze": "", "period": ""}],
                   out, "freeze_events.csv")
        save_table([{k: (v.isoformat() if isinstance(v, date) else v) for k, v in f.items()} for f in filled]
                   or [{"date": "", "variable": "", "value_c": "", "gap_days": ""}], out, "filled_days.csv")
        figures = make_figures(daily, stages, base_c, upper_c if method == "modified" else None,
                               args.freeze_c, out)
        ctx = {
            "warnings": warnings, "lat": args.lat, "lon": args.lon, "source_desc": source_desc,
            "crop": args.crop, "plant": plant.isoformat(), "end": end_eff.isoformat(),
            "n_days": len(daily), "method": method, "base_c": base_c, "upper_c": upper_c,
            "cum_c": daily[-1]["cum_gdd_c"], "cum_f": daily[-1]["cum_gdd_f"],
            "mean_tmax": sum(d["tmax_c"] for d in daily) / len(daily),
            "mean_tmin": sum(d["tmin_c"] for d in daily) / len(daily),
            "n_filled_days": len({f["date"] for f in filled}), "stages": stages,
            "hybrid": args.hybrid_gdd, "frost": frost, "freeze_c": args.freeze_c,
            "hard_c": args.hard_freeze_c, "figures": figures,
        }
        write_report(out, ctx)
        write_repro_bundle(out, argv, inputs=inputs, packages=["matplotlib"])
    except (ValueError, RuntimeError) as e:
        print(ascii_safe(f"Error: {e}"), file=sys.stderr)
        return 1
    msg = (f"Done: {len(daily)} days, {daily[-1]['cum_gdd_f']:.0f} degF-days "
           f"({daily[-1]['cum_gdd_c']:.0f} degC-days).")
    if stages:
        r1 = next(s for s in stages if s["stage"] == "R1")
        msg += f" R1: {r1['date_reached'] or 'not reached'}."
    print(ascii_safe(msg + f" Report: {(Path(args.out) / 'report.md').resolve()}"))
    for w in warnings:
        print(ascii_safe(f"Warning: {w}"), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
