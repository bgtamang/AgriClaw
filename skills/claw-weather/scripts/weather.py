#!/usr/bin/env python3
"""claw-weather: growing-season weather summary (temperature, rain, reference ET0,
simple climatic water balance, sourced heat/frost/heavy-rain indices, dry spells)
from Open-Meteo or a station CSV, with missing data counted, never treated as zero.

Usage:
    python scripts/weather.py --demo
    python scripts/weather.py --lat 40.08 --lon -88.24 --start-date 2025-04-01 --end-date 2025-10-31
    python scripts/weather.py --input station.csv --temp-units C --precip-units mm --start-date 2025-04-01

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

SKILL = "claw-weather"
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

ARCHIVE_URL = "https://archive-api.open-meteo.com/v1/archive"
FORECAST_URL = "https://api.open-meteo.com/v1/forecast"
PAST_DAYS_MAX = 92
MODELS = ("best_match", "era5", "era5_seamless", "era5_land", "ecmwf_ifs")
# Open-Meteo daily variable -> (output column, unit expected in daily_units)
OM_VARS = {
    "temperature_2m_max": ("tmax_c", "°C"),
    "temperature_2m_min": ("tmin_c", "°C"),
    "precipitation_sum": ("precip_mm", "mm"),
    "et0_fao_evapotranspiration": ("et0_mm", "mm"),
    "relative_humidity_2m_min": ("rh_min_pct", "%"),
    "wind_speed_10m_mean": ("wind_10m_m_s", "m/s"),
    "shortwave_radiation_sum": ("srad_mj_m2", "MJ/m²"),
}
DAILY_COLS = ["tmax_c", "tmin_c", "precip_mm", "et0_mm", "rh_min_pct", "wind_2m_m_s", "srad_mj_m2"]

WET_DAY_MM = 1.0                 # ETCCDI wet day: RR >= 1 mm (dry day: RR < 1 mm)
HARD_FREEZE_C = (28.0 - 32.0) / 1.8
HOT_C = (90.0 - 32.0) / 1.8      # NCH-40: rain-fed corn begins to stress above 90 degF
VERY_HOT_C = 35.0                # NCH-40: yield loss for each day at or above 95 degF (35 degC)
MIN_DRY_SPELL = 7                # listing threshold only (adjustable rule of thumb)

DEMO = {"lat": 40.08, "lon": -88.24, "start": "2025-04-01", "end": "2025-10-31",
        "file": "openmeteo_champaign_2025.json"}


# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------


def parse_date(text: str, what: str = "date") -> date:
    try:
        return date.fromisoformat(str(text).strip())
    except ValueError:
        raise ValueError(f"{what} '{text}' is not a valid YYYY-MM-DD date.") from None


def format_latlon(lat: float | None, lon: float | None) -> str:
    """40.08, -88.24 -> '40.0800°N, 88.2400°W'."""
    if lat is None or lon is None:
        return "not given"
    return (f"{abs(lat):.4f}°{'N' if lat >= 0 else 'S'}, "
            f"{abs(lon):.4f}°{'E' if lon >= 0 else 'W'}")


def wind_to_2m(u_z: float, z: float = 10.0) -> float:
    """FAO-56 eq. 47: u2 = uz * 4.87 / ln(67.8 z - 5.42) (log profile over short grass)."""
    return u_z * 4.87 / math.log(67.8 * z - 5.42)


def c_to_f(t: float) -> float:
    return t * 1.8 + 32.0


def f_to_c(t: float) -> float:
    return (t - 32.0) / 1.8


# ---------------------------------------------------------------------------
# Open-Meteo access (fetching is separate from parsing so tests run offline)
# ---------------------------------------------------------------------------


def openmeteo_url(kind: str, lat: float, lon: float, start: date | None = None,
                  end: date | None = None, model: str = "best_match",
                  past_days: int | None = None) -> str:
    params = {"latitude": f"{lat:.5f}", "longitude": f"{lon:.5f}", "daily": ",".join(OM_VARS),
              "timezone": "auto", "wind_speed_unit": "ms"}
    if kind == "archive":
        params.update(start_date=start.isoformat(), end_date=end.isoformat())
        if model != "best_match":
            params["models"] = model
        return ARCHIVE_URL + "?" + urllib.parse.urlencode(params)
    params.update(past_days=str(past_days), forecast_days="0")
    return FORECAST_URL + "?" + urllib.parse.urlencode(params)


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
        raise RuntimeError(f"Could not reach Open-Meteo ({e}). Check the connection, or use "
                           "--input with station data.") from None


def parse_openmeteo_daily(resp: dict, source: str) -> tuple[list[dict], dict]:
    if resp.get("error"):
        raise RuntimeError(f"Open-Meteo error: {resp.get('reason', resp)}")
    daily = resp.get("daily")
    if not isinstance(daily, dict) or "time" not in daily:
        raise RuntimeError("Open-Meteo response has no 'daily' block; cannot parse it.")
    units = resp.get("daily_units", {})
    rows = []
    for var, (_, unit) in OM_VARS.items():
        if var in daily and units.get(var) not in (None, unit):
            raise RuntimeError(f"'{var}' is in {units.get(var)}, expected {unit}.")
    for i, t in enumerate(daily["time"]):
        row = {"date": parse_date(t)}
        for var, (col, _) in OM_VARS.items():
            v = daily.get(var, [None] * len(daily["time"]))[i]
            row[col] = None if v is None else float(v)
        w = row.pop("wind_10m_m_s")
        row["wind_2m_m_s"] = None if w is None else wind_to_2m(w, 10.0)
        row["source"] = source
        rows.append(row)
    meta = {k: resp.get(k) for k in ("latitude", "longitude", "elevation", "timezone",
                                     "utc_offset_seconds")}
    return rows, meta


def local_today(utc_offset_seconds: int | None, now_utc: datetime | None = None) -> date:
    now_utc = now_utc or datetime.now(timezone.utc)
    return (now_utc + timedelta(seconds=utc_offset_seconds or 0)).date()


def trim_trailing_missing(rows: list[dict]) -> tuple[list[dict], int]:
    n = len(rows)
    while n and all(rows[n - 1][c] is None for c in DAILY_COLS):
        n -= 1
    return rows[:n], len(rows) - n


def fetch_daily(lat: float, lon: float, start: date, end: date, model: str = "best_match",
                today: date | None = None, fetcher=http_get_json,
                now_utc: datetime | None = None) -> tuple[list[dict], dict, list[str], list[dict]]:
    """Archive data for [start, end] with the end clamped to yesterday; days the
    archive has not published yet are filled from the forecast API (&past_days=)."""
    today = today or date.today()
    notes, raw = [], []
    yesterday = today - timedelta(days=1)
    if start > yesterday:
        raise ValueError(f"Start date {start} is not in the past (today is {today}).")
    end_eff = min(end, yesterday)
    if end_eff < end:
        notes.append(f"Requested end date {end} is after the last complete day; the summary "
                     f"stops at {end_eff} at the latest. Nothing is forecast.")
    url = openmeteo_url("archive", lat, lon, start, end_eff, model=model)
    resp = fetcher(url)
    raw.append({"url": url, "response": resp})
    rows, meta = parse_openmeteo_daily(resp, f"archive:{model}")
    ltoday = local_today(meta.get("utc_offset_seconds"), now_utc)
    rows = [r for r in rows if r["date"] < ltoday]
    rows, n_trim = trim_trailing_missing(rows)
    if n_trim:
        notes.append(f"The archive ({model}) has no data yet for the last {n_trim} day(s) requested "
                     "(reanalysis is published with a delay).")
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
                notes.append(f"{len(rows2)} most recent day(s) ({rows2[0]['date']} to {rows2[-1]['date']}) "
                             "come from the Open-Meteo forecast API (&past_days=): operational model "
                             "output, provisional, not the archive dataset.")
        else:
            notes.append(f"No data after {last}: not yet in the archive and older than the forecast "
                         "API's 92-day past_days window.")
    return rows, meta, notes, raw


# ---------------------------------------------------------------------------
# Station CSV input
# ---------------------------------------------------------------------------

MISSING_TOKENS = {"", "na", "nan", "null", "none", "missing"}


def read_station_csv(path: str, temp_units: str | None, precip_units: str | None) -> list[dict]:
    """date plus any of tmax, tmin (--temp-units), precip, et0 (--precip-units),
    rh_min_pct, wind_2m_m_s or wind_10m_m_s, srad_mj_m2. Units are never inferred."""
    p = Path(path)
    if not p.is_file():
        raise ValueError(f"Input file not found: {path}")
    with p.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        header = {h.strip().lower(): h for h in (reader.fieldnames or [])}
        raw = list(reader)
    if "date" not in header:
        raise ValueError(f"No 'date' column in {path}; found {', '.join(header) or 'none'}.")
    spec = {}  # output column -> (input header, converter)
    for key, units, opts in (("tmax", temp_units, ("c", "f")), ("tmin", temp_units, ("c", "f")),
                             ("precip", precip_units, ("mm", "in")), ("et0", precip_units, ("mm", "in"))):
        found = [h for h in [key] + [f"{key}_{o}" for o in opts] if h in header]
        if not found:
            continue
        if len(found) > 1:
            raise ValueError(f"Columns {found} are ambiguous; keep one {key} column.")
        h = found[0]
        flag = "--temp-units" if key in ("tmax", "tmin") else "--precip-units"
        if units is None:
            raise ValueError(f"Column '{header[h]}' needs {flag} (units are never guessed).")
        suffix = h.split("_")[-1] if "_" in h else None
        if suffix and suffix != units.lower():
            raise ValueError(f"Column '{header[h]}' says {suffix} but {flag} is {units}.")
        if key in ("tmax", "tmin"):
            conv = f_to_c if units == "F" else (lambda v: v)
            spec[f"{key}_c"] = (header[h], conv)
        else:
            conv = (lambda v: v * 25.4) if units == "in" else (lambda v: v)
            spec[f"{key}_mm"] = (header[h], conv)
    if "rh_min_pct" in header:
        spec["rh_min_pct"] = (header["rh_min_pct"], lambda v: v)
    if "wind_2m_m_s" in header:
        spec["wind_2m_m_s"] = (header["wind_2m_m_s"], lambda v: v)
    elif "wind_10m_m_s" in header:
        spec["wind_2m_m_s"] = (header["wind_10m_m_s"], lambda v: wind_to_2m(v, 10.0))
    if "srad_mj_m2" in header:
        spec["srad_mj_m2"] = (header["srad_mj_m2"], lambda v: v)
    if not spec:
        raise ValueError("The CSV has no recognized weather columns (tmax, tmin, precip, et0, ...).")
    rows, seen = [], set()
    for i, r in enumerate(raw, start=2):
        d = parse_date(r[header["date"]], f"Row {i}: date")
        if d in seen:
            raise ValueError(f"Row {i}: date {d} appears more than once.")
        seen.add(d)
        row = {"date": d}
        for col in DAILY_COLS:
            if col not in spec:
                row[col] = None
                continue
            h, conv = spec[col]
            txt = (r[h] or "").strip()
            if txt.lower() in MISSING_TOKENS:
                row[col] = None
                continue
            try:
                row[col] = conv(float(txt))
            except ValueError:
                raise ValueError(f"Row {i}: {h} '{txt}' is not numeric.") from None
        row["source"] = "input CSV"
        rows.append(row)
    rows.sort(key=lambda r: r["date"])
    check_sanity(rows)
    return rows


def check_sanity(rows: list[dict]) -> None:
    checks = [("tmax_c", -60, 60), ("tmin_c", -60, 60), ("precip_mm", 0, 1000), ("et0_mm", 0, 20),
              ("rh_min_pct", 0, 100), ("wind_2m_m_s", 0, 60), ("srad_mj_m2", 0, 50)]
    for col, lo, hi in checks:
        bad = [r for r in rows if r[col] is not None and not lo <= r[col] <= hi]
        if bad:
            raise ValueError(f"{len(bad)} value(s) of {col} outside {lo}..{hi} (first {bad[0]['date']}: "
                             f"{bad[0][col]}). Check the units flags and replace missing-value codes "
                             "such as -9999 with blank cells.")


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------


def complete_calendar(rows: list[dict], start: date, end: date) -> list[dict]:
    by_date = {r["date"]: r for r in rows}
    out, d = [], start
    while d <= end:
        out.append(dict(by_date.get(d) or dict({c: None for c in DAILY_COLS}, date=d, source="missing")))
        d += timedelta(days=1)
    return out


def valid(rows: list[dict], col: str) -> list[float]:
    return [r[col] for r in rows if r[col] is not None]


def summarize(rows: list[dict]) -> dict:
    """Totals and means over days with data; missing days are counted, not zeroed."""
    s = {"n_days": len(rows)}
    for col in DAILY_COLS:
        v = valid(rows, col)
        s[f"n_missing_{col}"] = len(rows) - len(v)
    for col in ("precip_mm", "et0_mm", "srad_mj_m2"):
        v = valid(rows, col)
        s[f"total_{col}"] = round(sum(v), 1) if v else None
    for col in ("tmax_c", "tmin_c", "rh_min_pct", "wind_2m_m_s"):
        v = valid(rows, col)
        s[f"mean_{col}"] = round(sum(v) / len(v), 2) if v else None
    both = [(r["precip_mm"], r["et0_mm"]) for r in rows
            if r["precip_mm"] is not None and r["et0_mm"] is not None]
    s["p_minus_et0_mm"] = round(sum(p - e for p, e in both), 1) if both else None
    s["n_days_p_and_et0"] = len(both)
    tx = [(r["tmax_c"], r["date"]) for r in rows if r["tmax_c"] is not None]
    tn = [(r["tmin_c"], r["date"]) for r in rows if r["tmin_c"] is not None]
    s["max_tmax_c"], s["date_max_tmax"] = (round(max(tx)[0], 1), max(tx)[1].isoformat()) if tx else (None, "")
    s["min_tmin_c"], s["date_min_tmin"] = (round(min(tn)[0], 1), min(tn)[1].isoformat()) if tn else (None, "")
    core = ("tmax_c", "tmin_c", "precip_mm", "et0_mm")
    s["complete"] = "yes" if all(s[f"n_missing_{c}"] == 0 for c in core) else "no"
    return s


def monthly(rows: list[dict]) -> list[dict]:
    months: dict[str, list[dict]] = {}
    for r in rows:
        months.setdefault(r["date"].strftime("%Y-%m"), []).append(r)
    out = []
    for m in sorted(months):
        s = summarize(months[m])
        out.append({"month": m, "n_days": s["n_days"],
                    "mean_tmax_c": s["mean_tmax_c"], "mean_tmin_c": s["mean_tmin_c"],
                    "precip_mm": s["total_precip_mm"], "et0_mm": s["total_et0_mm"],
                    "p_minus_et0_mm": s["p_minus_et0_mm"], "srad_mj_m2": s["total_srad_mj_m2"],
                    "n_missing_tmax": s["n_missing_tmax_c"], "n_missing_tmin": s["n_missing_tmin_c"],
                    "n_missing_precip": s["n_missing_precip_mm"], "n_missing_et0": s["n_missing_et0_mm"],
                    "complete": s["complete"]})
    return out


def threshold_indices(rows: list[dict], hot_c: float = HOT_C, very_hot_c: float = VERY_HOT_C,
                      hard_freeze_c: float = HARD_FREEZE_C, extra_rain_mm: float | None = None):
    """Counts over the analyzed period. Each index uses one variable; days where that
    variable is missing are counted separately (n_missing)."""
    defs = [
        ("FD", "frost days", "tmin_c", lambda v: v < 0.0, "TN < 0 °C", "ETCCDI (Zhang et al. 2011)"),
        ("HF28", "hard-freeze days", "tmin_c", lambda v: v <= hard_freeze_c,
         f"TN ≤ {hard_freeze_c:.2f} °C (28 °F)", "NOAA NCEI hard freeze (Arguez 2012)"),
        ("SU", "summer days", "tmax_c", lambda v: v > 25.0, "TX > 25 °C", "ETCCDI (Zhang et al. 2011)"),
        ("TR", "tropical nights", "tmin_c", lambda v: v > 20.0, "TN > 20 °C", "ETCCDI (Zhang et al. 2011)"),
        ("TX90F", "hot days (corn stress)", "tmax_c", lambda v: v > hot_c,
         f"TX > {hot_c:.1f} °C ({c_to_f(hot_c):.0f} °F)", "Neild & Newman (1987) NCH-40"),
        ("TX95F", "very hot days (corn yield loss)", "tmax_c", lambda v: v >= very_hot_c,
         f"TX ≥ {very_hot_c:.1f} °C ({c_to_f(very_hot_c):.0f} °F)", "Neild & Newman (1987) NCH-40"),
        ("R10mm", "heavy-precipitation days", "precip_mm", lambda v: v >= 10.0, "RR ≥ 10 mm",
         "ETCCDI (Zhang et al. 2011)"),
        ("R20mm", "very-heavy-precipitation days", "precip_mm", lambda v: v >= 20.0, "RR ≥ 20 mm",
         "ETCCDI (Zhang et al. 2011)"),
    ]
    if extra_rain_mm is not None:
        defs.append((f"R{extra_rain_mm:g}mm", "user-defined heavy-rain days", "precip_mm",
                     lambda v, t=extra_rain_mm: v >= t, f"RR ≥ {extra_rain_mm:g} mm",
                     "ETCCDI Rnnmm with a user threshold"))
    table, events = [], []
    for code, name, col, test, rule, source in defs:
        hits = [r for r in rows if r[col] is not None and test(r[col])]
        table.append({"index": code, "name": name, "rule": rule, "days": len(hits),
                      "n_missing": sum(r[col] is None for r in rows), "source": source})
        events += [{"date": r["date"].isoformat(), "index": code, "variable": col,
                    "value": round(r[col], 2)} for r in hits]
    p = [r for r in rows if r["precip_mm"] is not None]
    rx1 = max(p, key=lambda r: r["precip_mm"]) if p else None
    table.append({"index": "Rx1day", "name": "maximum 1-day precipitation", "rule": "max RR (mm)",
                  "days": round(rx1["precip_mm"], 1) if rx1 else "",
                  "n_missing": sum(r["precip_mm"] is None for r in rows), "source": "ETCCDI (Zhang et al. 2011)"})
    spells = dry_spells(rows)
    wet = wet_spells(rows)
    table.append({"index": "CDD", "name": "longest dry spell (days)", "rule": "max run of RR < 1 mm",
                  "days": max((s["length_days"] for s in spells), default=0),
                  "n_missing": sum(r["precip_mm"] is None for r in rows), "source": "ETCCDI (Zhang et al. 2011)"})
    table.append({"index": "CWD", "name": "longest wet spell (days)", "rule": "max run of RR ≥ 1 mm",
                  "days": max(wet, default=0),
                  "n_missing": sum(r["precip_mm"] is None for r in rows), "source": "ETCCDI (Zhang et al. 2011)"})
    events.sort(key=lambda e: (e["date"], e["index"]))
    return table, events, spells


def dry_spells(rows: list[dict], wet_mm: float = WET_DAY_MM) -> list[dict]:
    """Runs of consecutive days with RR < 1 mm. A missing day ends a run; runs that
    touch missing data or the period edges are flagged (true length may be longer)."""
    spells, start, n = [], None, 0

    def close(i_end):
        before = rows[i_end - n - 1]["precip_mm"] if i_end - n - 1 >= 0 else "edge"
        after = rows[i_end]["precip_mm"] if i_end < len(rows) else "edge"
        spells.append({"start": rows[i_end - n]["date"].isoformat(), "end": rows[i_end - 1]["date"].isoformat(),
                       "length_days": n,
                       "touches_missing": "yes" if before is None or after is None else "no",
                       "touches_period_edge": "yes" if "edge" in (before, after) else "no"})

    for i, r in enumerate(rows):
        v = r["precip_mm"]
        if v is not None and v < wet_mm:
            n += 1
        else:
            if n:
                close(i)
            n = 0
    if n:
        close(len(rows))
    spells.sort(key=lambda s: (-s["length_days"], s["start"]))
    return spells


def wet_spells(rows: list[dict], wet_mm: float = WET_DAY_MM) -> list[int]:
    out, n = [], 0
    for r in rows:
        v = r["precip_mm"]
        if v is not None and v >= wet_mm:
            n += 1
        else:
            if n:
                out.append(n)
            n = 0
    if n:
        out.append(n)
    return out


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def make_figures(rows: list[dict], months: list[dict], out: Path) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    fig_dir = out / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    dates = [r["date"] for r in rows]
    nan = float("nan")
    files = []

    fig, ax1 = plt.subplots(figsize=(10, 4.5))
    ax2 = ax1.twinx()
    ax2.bar(dates, [r["precip_mm"] if r["precip_mm"] is not None else 0 for r in rows],
            width=1.0, color="#4a90c2", alpha=0.6, label="precipitation (mm)")
    miss = [r["date"] for r in rows if r["precip_mm"] is None]
    if miss:
        ax2.plot(miss, [0] * len(miss), "|", color="red", ms=10, label="precipitation missing")
    ax2.set_ylabel("Precipitation (mm/day)")
    ax1.plot(dates, [r["tmax_c"] if r["tmax_c"] is not None else nan for r in rows], color="#c0392b", lw=1,
             label="Tmax")
    ax1.plot(dates, [r["tmin_c"] if r["tmin_c"] is not None else nan for r in rows], color="#2e6da4", lw=1,
             label="Tmin")
    ax1.set_ylabel("Air temperature (°C)")
    ax1.set_zorder(ax2.get_zorder() + 1)
    ax1.patch.set_visible(False)
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper left")
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.set_title("Daily temperature and precipitation")
    ax1.grid(alpha=0.3)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(fig_dir / "temp_precip.png", dpi=150)
    plt.close(fig)
    files.append("figures/temp_precip.png")

    if any(r["et0_mm"] is not None for r in rows):
        cp, ce, cpe = [], [], []
        sp = se = 0.0
        for r in rows:
            sp += 0.0 if r["precip_mm"] is None else r["precip_mm"]
            se += 0.0 if r["et0_mm"] is None else r["et0_mm"]
            cp.append(sp)
            ce.append(se)
        fig, ax = plt.subplots(figsize=(10, 4.5))
        ax.plot(dates, cp, color="#2e6da4", lw=2, label="cumulative precipitation")
        ax.plot(dates, ce, color="#c0392b", lw=2, label="cumulative reference ET0")
        ax.fill_between(dates, cp, ce, where=[a < b for a, b in zip(cp, ce)], color="#c0392b", alpha=0.15,
                        label="ET0 > P")
        ax.set_ylabel("mm (missing days add nothing)")
        ax.set_title("Simple climatic water balance (not crop water use)")
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(alpha=0.3)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        fig.autofmt_xdate()
        fig.tight_layout()
        fig.savefig(fig_dir / "water_balance.png", dpi=150)
        plt.close(fig)
        files.append("figures/water_balance.png")

    fig, ax = plt.subplots(figsize=(9, 4.5))
    x = list(range(len(months)))
    ax.bar([i - 0.2 for i in x], [m["precip_mm"] or 0 for m in months], 0.4, color="#4a90c2", label="P (mm)")
    ax.bar([i + 0.2 for i in x], [m["et0_mm"] or 0 for m in months], 0.4, color="#e8a39b", label="ET0 (mm)")
    for i, m in enumerate(months):
        if m["complete"] == "no":
            ax.annotate("incomplete", (i, 0), xytext=(0, 3), textcoords="offset points", ha="center",
                        fontsize=7, color="red", rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels([m["month"] for m in months], rotation=45, ha="right")
    ax.set_ylabel("mm per month")
    ax2 = ax.twinx()
    ax2.plot(x, [m["mean_tmax_c"] if m["mean_tmax_c"] is not None else nan for m in months], "o-",
             color="#c0392b", label="mean Tmax")
    ax2.plot(x, [m["mean_tmin_c"] if m["mean_tmin_c"] is not None else nan for m in months], "s-",
             color="#2e6da4", label="mean Tmin")
    ax2.set_ylabel("°C")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper left")
    ax.set_title("Monthly summary")
    fig.tight_layout()
    fig.savefig(fig_dir / "monthly_summary.png", dpi=150)
    plt.close(fig)
    files.append("figures/monthly_summary.png")
    return files


REFERENCES = [
    "Zhang, X., Alexander, L., Hegerl, G.C., Jones, P., Klein Tank, A., Peterson, T.C., Trewin, B. & "
    "Zwiers, F.W. (2011). Indices for monitoring changes in extremes based on daily temperature and "
    "precipitation data. WIREs Climate Change 2(6), 851-870. https://doi.org/10.1002/wcc.147 "
    "(definitions: https://etccdi.pacificclimate.org/list_27_indices.shtml)",
    "Neild, R.E. & Newman, J.E. (1987; Purdue reprint RR 4/90). Growing Season Characteristics and "
    "Requirements in the Corn Belt. National Corn Handbook NCH-40. "
    "https://www.extension.purdue.edu/extmedia/nch/nch-40.html",
    "Arguez, A. (2012). NOAA's 1981-2010 Climate Normals: Supplemental Normals. NOAA NCDC technical report, "
    "in NOAA NCEI Normals Calculation Methodology 2020. https://www.ncei.noaa.gov/data/"
    "normals-annualseasonal/1991-2020/doc/Normals_Calculation_Methodology_2020.pdf",
    "Allen, R.G., Pereira, L.S., Raes, D. & Smith, M. (1998). Crop evapotranspiration: guidelines for "
    "computing crop water requirements. FAO Irrigation and Drainage Paper 56. FAO, Rome. "
    "https://www.fao.org/4/x0490e/x0490e00.htm",
    "Zippenfenig, P. (2023). Open-Meteo.com Weather API [Computer software]. Zenodo. "
    "https://doi.org/10.5281/zenodo.7970649",
    "Hersbach, H. et al. (2020). The ERA5 global reanalysis. Quarterly Journal of the Royal "
    "Meteorological Society 146(730), 1999-2049. https://doi.org/10.1002/qj.3803",
    "Muñoz-Sabater, J. et al. (2021). ERA5-Land: a state-of-the-art global reanalysis dataset for land "
    "applications. Earth System Science Data 13(9), 4349-4383. https://doi.org/10.5194/essd-13-4349-2021",
]


def fmt(v, nd=1, unit=""):
    if v is None or v == "":
        return "n/a"
    return f"{v:.{nd}f}{unit}" if isinstance(v, float) else f"{v}{unit}"


def write_report(out: Path, ctx: dict) -> None:
    s = ctx["summary"]
    L = [f"# {SKILL} report", ""]
    if ctx["warnings"]:
        L += ["## Warnings", ""] + [f"- {w}" for w in ctx["warnings"]] + [""]
    L += ["## Data source", "",
          f"- Location: {format_latlon(ctx['lat'], ctx['lon'])}",
          f"- Source: {ctx['source_desc']}",
          f"- Period: {ctx['start']} to {ctx['end']} ({s['n_days']} days)", ""]
    if ctx["gridded"]:
        L += ["Open-Meteo's archive is gridded model output, not station observations. The default 'best match' "
              "combines ECMWF IFS operational model data (9 km, from 2017) with the ERA5 (0.25°, ~25 km) and "
              "ERA5-Land (0.1°, ~11 km) reanalyses. One grid-cell value can differ from a field's own weather, especially for "
              "convective rain, frost in low spots and near coasts or mountains. ET0 is Open-Meteo's FAO-56 "
              "Penman-Monteith grass reference ET.", ""]
    L += ["## Season summary", "", "| Metric | Value | Missing days |", "|---|---|---|",
          f"| Total precipitation | {fmt(s['total_precip_mm'], 1, ' mm')} | {s['n_missing_precip_mm']} |",
          f"| Total reference ET0 | {fmt(s['total_et0_mm'], 1, ' mm')} | {s['n_missing_et0_mm']} |",
          f"| P − ET0 (days with both) | {fmt(s['p_minus_et0_mm'], 1, ' mm')} | {s['n_days'] - s['n_days_p_and_et0']} |",
          f"| Mean Tmax / Tmin | {fmt(s['mean_tmax_c'], 1)} / {fmt(s['mean_tmin_c'], 1)} °C | "
          f"{s['n_missing_tmax_c']} / {s['n_missing_tmin_c']} |",
          f"| Highest Tmax | {fmt(s['max_tmax_c'], 1, ' °C')} ({s['date_max_tmax']}) | |",
          f"| Lowest Tmin | {fmt(s['min_tmin_c'], 1, ' °C')} ({s['date_min_tmin']}) | |",
          f"| Total shortwave radiation | {fmt(s['total_srad_mj_m2'], 0, ' MJ/m²')} | {s['n_missing_srad_mj_m2']} |",
          "", "Totals add up the days that have data; a total with missing days is a lower bound, "
          "not an estimate of the full period.", ""]
    L += ["## Monthly summary and simple climatic water balance", "",
          "P − ET0 compares rainfall with the evaporative demand of a reference grass surface. It is a simple "
          "climatic water balance, **not a drought index and not crop water use**: it ignores soil water "
          "storage, runoff, drainage and the crop's growth stage. Use **claw-water-balance** for crop "
          "evapotranspiration and water stress.", "",
          "| Month | Days | Mean Tmax (°C) | Mean Tmin (°C) | P (mm) | ET0 (mm) | P − ET0 (mm) | Complete |",
          "|---|---|---|---|---|---|---|---|"]
    for m in ctx["months"]:
        L.append(f"| {m['month']} | {m['n_days']} | {fmt(m['mean_tmax_c'])} | {fmt(m['mean_tmin_c'])} | "
                 f"{fmt(m['precip_mm'])} | {fmt(m['et0_mm'])} | {fmt(m['p_minus_et0_mm'])} | {m['complete']} |")
    L += ["", "## Threshold indices", "",
          "Counts over the analyzed period (ETCCDI defines them per year; here they cover only "
          f"{ctx['start']} to {ctx['end']}). Days with the variable missing are not counted.", "",
          "| Index | Meaning | Rule | Days | Missing days | Source |", "|---|---|---|---|---|---|"]
    for t in ctx["indices"]:
        L.append(f"| {t['index']} | {t['name']} | {t['rule']} | {t['days']} | {t['n_missing']} | {t['source']} |")
    L += ["", f"## Dry spells (runs of days with < {WET_DAY_MM:g} mm, ≥ {ctx['min_spell']} days)", ""]
    longs = [sp for sp in ctx["spells"] if sp["length_days"] >= ctx["min_spell"]]
    if longs:
        L += ["| Start | End | Days | Touches missing data | Touches period edge |", "|---|---|---|---|---|"]
        L += [f"| {sp['start']} | {sp['end']} | {sp['length_days']} | {sp['touches_missing']} | "
              f"{sp['touches_period_edge']} |" for sp in longs[:15]]
    else:
        L.append(f"No dry spell of {ctx['min_spell']} days or more.")
    L += ["", "## Figures", ""] + [f"![{f}]({f})" for f in ctx["figures"]] + [""]
    L += ["## References", ""] + [f"- {r}" for r in REFERENCES] + ["", f"*{SKILL} {VERSION}*"]
    write_text(out / "report.md", "\n".join(L) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="weather.py", description=__doc__.splitlines()[0])
    ap.add_argument("--demo", action="store_true",
                    help="Champaign, IL, April 1 - October 31, 2025 (bundled Open-Meteo data)")
    ap.add_argument("--lat", type=float, help="latitude, decimal degrees")
    ap.add_argument("--lon", type=float, help="longitude, decimal degrees")
    ap.add_argument("--start-date", help="first day, YYYY-MM-DD (required unless --demo)")
    ap.add_argument("--end-date", help="last day, YYYY-MM-DD (default: yesterday, or the last CSV date)")
    ap.add_argument("--input", help="station CSV (date plus tmax/tmin/precip/et0/...)")
    ap.add_argument("--temp-units", choices=["C", "F"], help="units of tmax/tmin in --input")
    ap.add_argument("--precip-units", choices=["mm", "in"], help="units of precip/et0 in --input")
    ap.add_argument("--model", choices=MODELS, default="best_match", help="Open-Meteo archive model")
    ap.add_argument("--hot-c", type=float, default=HOT_C, help="hot-day threshold, TX > value (°C)")
    ap.add_argument("--very-hot-c", type=float, default=VERY_HOT_C, help="very-hot threshold, TX ≥ value (°C)")
    ap.add_argument("--heavy-rain-mm", type=float, help="extra user-defined heavy-rain threshold (mm)")
    ap.add_argument("--min-dry-spell", type=int, default=MIN_DRY_SPELL, help="list dry spells of at least N days")
    ap.add_argument("--out", default=f"agriclaw_output/{SKILL}", help="output folder")
    return ap


def main(argv: list[str] | None = None, today: date | None = None, fetcher=http_get_json) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = build_parser()
    args = ap.parse_args(argv)
    warnings: list[str] = []
    try:
        if args.demo:
            args.lat = DEMO["lat"] if args.lat is None else args.lat
            args.lon = DEMO["lon"] if args.lon is None else args.lon
            args.start_date = args.start_date or DEMO["start"]
            args.end_date = args.end_date or DEMO["end"]
        if not args.start_date:
            ap.error("--start-date is required (ask the user for the period to summarize)")
        if not (args.demo or args.input) and (args.lat is None or args.lon is None):
            ap.error("give --lat and --lon, or --input FILE.csv, or --demo")
        for v, lo, hi, name in ((args.lat, -90, 90, "--lat"), (args.lon, -180, 180, "--lon")):
            if v is not None and not lo <= v <= hi:
                ap.error(f"{name} must be within {lo}..{hi}")
        start = parse_date(args.start_date, "--start-date")
        today = today or date.today()
        out = Path(args.out)
        out.mkdir(parents=True, exist_ok=True)
        inputs, gridded = [], True
        if args.input:
            rows = read_station_csv(args.input, args.temp_units, args.precip_units)
            if not rows:
                raise ValueError(f"{args.input} has no data rows.")
            inputs.append(args.input)
            gridded = False
            end = parse_date(args.end_date, "--end-date") if args.end_date else rows[-1]["date"]
            source_desc = f"station CSV `{args.input}`"
            if rows[0]["date"] > start or rows[-1]["date"] < end:
                warnings.append(f"The CSV covers {rows[0]['date']} to {rows[-1]['date']}; days of the requested "
                                f"period outside that range are counted as missing.")
        elif args.demo:
            f = SKILL_DIR / "demo" / DEMO["file"]
            rows, meta = parse_openmeteo_daily(json.loads(f.read_text(encoding="utf-8")),
                                               "archive:best_match (bundled demo file)")
            inputs.append(str(f))
            end = parse_date(args.end_date, "--end-date")
            source_desc = (f"Open-Meteo archive, best match, grid cell {format_latlon(meta['latitude'], meta['longitude'])}"
                           f" (saved in demo/{DEMO['file']})")
        else:
            end = parse_date(args.end_date, "--end-date") if args.end_date else today - timedelta(days=1)
            if end < start:
                raise ValueError("--end-date is before --start-date.")
            rows, meta, notes, raw = fetch_daily(args.lat, args.lon, start, end, args.model, today, fetcher)
            warnings += notes
            write_text(out / "raw" / "openmeteo_responses.json", json.dumps(raw, indent=1))
            end = min(end, rows[-1]["date"]) if rows else start
            source_desc = (f"Open-Meteo archive, model '{args.model}', grid cell "
                           f"{format_latlon(meta['latitude'], meta['longitude'])}")
        if end < start:
            raise ValueError("--end-date is before --start-date.")
        rows = complete_calendar([r for r in rows if start <= r["date"] <= end], start, end)
        if not any(r[c] is not None for r in rows for c in DAILY_COLS):
            raise ValueError("No weather data in the requested period.")
        for col in DAILY_COLS:
            n_miss = sum(r[col] is None for r in rows)
            if n_miss == len(rows):
                warnings.append(f"{col} is not available from this source; it is left out of the summary.")
            elif n_miss:
                warnings.append(f"{col}: {n_miss} of {len(rows)} days missing. Totals and counts use only days "
                                "with data; missing days were never treated as zero.")
        if any(r["source"] == "forecast:past_days" for r in rows) and not any("forecast API" in w for w in warnings):
            warnings.append("Some days come from the forecast API (provisional).")

        summary = summarize(rows)
        months = monthly(rows)
        indices, events, spells = threshold_indices(rows, args.hot_c, args.very_hot_c, HARD_FREEZE_C,
                                                    args.heavy_rain_mm)
        if spells and spells[0]["touches_missing"] == "yes":
            warnings.append("The longest dry spell borders missing precipitation data; it may be longer.")

        daily_out = [{"date": r["date"].isoformat(),
                      **{c: (round(r[c], 2) if r[c] is not None else "") for c in DAILY_COLS},
                      "source": r["source"]} for r in rows]
        save_table(daily_out, out, "daily_weather.csv")
        save_table(months, out, "monthly_summary.csv")
        save_table([{k: v for k, v in summary.items()}], out, "season_summary.csv")
        save_table(indices, out, "threshold_indices.csv")
        save_table(events or [{"date": "", "index": "", "variable": "", "value": ""}], out, "threshold_days.csv")
        save_table(spells or [{"start": "", "end": "", "length_days": "", "touches_missing": "",
                               "touches_period_edge": ""}], out, "dry_spells.csv")
        partly = [c for c in DAILY_COLS if 0 < sum(r[c] is None for r in rows) < len(rows)]
        missing = [{"date": r["date"].isoformat(), "variable": c} for r in rows for c in partly if r[c] is None]
        save_table(missing or [{"date": "", "variable": ""}], out, "missing_days.csv")
        figures = make_figures(rows, months, out)
        write_report(out, {"warnings": warnings, "lat": args.lat, "lon": args.lon, "source_desc": source_desc,
                           "start": start.isoformat(), "end": end.isoformat(), "summary": summary,
                           "months": months, "indices": indices, "spells": spells,
                           "min_spell": args.min_dry_spell, "figures": figures, "gridded": gridded})
        write_repro_bundle(out, argv, inputs=inputs, packages=["matplotlib"])
    except (ValueError, RuntimeError) as e:
        print(ascii_safe(f"Error: {e}"), file=sys.stderr)
        return 1
    print(ascii_safe(f"Done: {summary['n_days']} days, precipitation {fmt(summary['total_precip_mm'])} mm, "
                     f"ET0 {fmt(summary['total_et0_mm'])} mm. Report: {(out / 'report.md').resolve()}"))
    for w in warnings:
        print(ascii_safe(f"Warning: {w}"), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
