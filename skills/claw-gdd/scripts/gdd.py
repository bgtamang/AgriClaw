#!/usr/bin/env python3
"""
claw-gdd: Growing Degree Day Calculator

Computes daily and cumulative GDD from weather data, estimates crop growth
stages, detects frost dates, and generates a reproducibility bundle.

Usage:
    python gdd.py --demo
    python gdd.py --lat 40.12 --lon -88.24 --plant-date 2025-05-15 --crop soybean
    python gdd.py --lat 40.12 --lon -88.24 --plant-date 2025-05-15 --crop soybean --base 8 --ceiling 35
    python gdd.py --input weather.csv --crop corn --plant-date 2025-04-20
"""

import argparse
import csv
import hashlib
import json
import os
import sys
from datetime import datetime, timedelta
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import URLError

# ---------------------------------------------------------------------------
# Crop parameters: base temp (C), ceiling temp (C), growth stage GDD thresholds
#
# References for default base and ceiling temperatures:
#
# Soybean (Tbase=10C, Tceiling=30C):
#   - Setiyono et al. (2007). Understanding and modeling the effect of temperature
#     and daylength on soybean phenology under high-yield conditions.
#     Field Crops Research, 100(2-3), 257-271.
#   - Gilmore & Rogers (1958). Heat units as a method of measuring maturity in corn.
#     Agronomy Journal, 50(10), 611-615. [foundational GDD method]
#   - University of Illinois Extension, "Growing Degree Days and Soybean Development"
#
# Corn (Tbase=10C, Tceiling=30C):
#   - Cross & Zuber (1972). Prediction of flowering dates in maize based on different
#     methods of estimating thermal units. Agronomy Journal, 64(3), 351-355.
#   - McMaster & Wilhelm (1997). Growing degree-days: one equation, two interpretations.
#     Agricultural and Forest Meteorology, 87(4), 291-300.
#   - USDA-NASS standard: base 50F (10C), ceiling 86F (30C)
#
# Wheat (Tbase=0C, Tceiling=25C):
#   - Porter & Gawith (1999). Temperatures and the growth and development of wheat:
#     a review. European Journal of Agronomy, 10(1), 23-36.
#   - McMaster & Smika (1988). Estimation and evaluation of winter wheat phenology
#     in the central Great Plains. Agricultural and Forest Meteorology, 43, 1-18.
#
# Rice (Tbase=10C, Tceiling=35C):
#   - Yoshida (1981). Fundamentals of Rice Crop Science. IRRI, Los Banos.
#   - Gao et al. (1992). An analysis of the base temperature for computing rice
#     growing degree days. Acta Agronomica Sinica, 18(1), 11-18.
#
# Growth stage GDD thresholds are approximate and vary by cultivar and maturity
# group. Values here represent mid-maturity cultivars in US Midwest conditions.
# Users should adjust with --base and --ceiling for their specific conditions.
# ---------------------------------------------------------------------------
CROP_PARAMS = {
    "soybean": {
        "base": 10.0,
        "ceiling": 30.0,
        "stages": [
            (90, "VE", "Emergence"),
            (130, "V1", "First trifoliate"),
            (280, "V3", "Third trifoliate"),
            (530, "V6", "Sixth trifoliate"),
            (800, "R1", "Beginning bloom"),
            (1050, "R3", "Beginning pod"),
            (1350, "R5", "Beginning seed fill"),
            (1650, "R7", "Beginning maturity"),
            (1800, "R8", "Full maturity"),
        ],
    },
    "corn": {
        "base": 10.0,
        "ceiling": 30.0,
        "stages": [
            (90, "VE", "Emergence"),
            (200, "V2", "Second leaf collar"),
            (475, "V6", "Sixth leaf collar"),
            (740, "V10", "Tenth leaf collar"),
            (1135, "VT", "Tasseling"),
            (1250, "R1", "Silking"),
            (1500, "R3", "Milk"),
            (1925, "R5", "Dent"),
            (2450, "R6", "Physiological maturity"),
        ],
    },
    "wheat": {
        "base": 0.0,
        "ceiling": 25.0,
        "stages": [
            (100, "Emergence", "Seedling emergence"),
            (400, "Tillering", "Tiller formation"),
            (800, "Jointing", "Stem elongation"),
            (1050, "Booting", "Flag leaf visible"),
            (1200, "Heading", "Head emergence"),
            (1400, "Flowering", "Anthesis"),
            (1700, "Grain fill", "Kernel development"),
            (2000, "Maturity", "Physiological maturity"),
        ],
    },
    "rice": {
        "base": 10.0,
        "ceiling": 35.0,
        "stages": [
            (100, "Emergence", "Seedling emergence"),
            (400, "Tillering", "Active tillering"),
            (800, "Panicle init.", "Panicle initiation"),
            (1100, "Booting", "Flag leaf visible"),
            (1300, "Heading", "Panicle emergence"),
            (1500, "Flowering", "Anthesis"),
            (1800, "Grain fill", "Grain development"),
            (2100, "Maturity", "Physiological maturity"),
        ],
    },
}

# ---------------------------------------------------------------------------
# Weather data fetching (Open-Meteo, free, no API key)
# ---------------------------------------------------------------------------

def fetch_weather(lat, lon, start_date, end_date):
    """Fetch daily Tmax/Tmin/precip from Open-Meteo historical archive."""
    url = (
        f"https://archive-api.open-meteo.com/v1/archive?"
        f"latitude={lat}&longitude={lon}"
        f"&start_date={start_date}&end_date={end_date}"
        f"&daily=temperature_2m_max,temperature_2m_min,precipitation_sum"
        f"&timezone=auto"
    )
    try:
        req = Request(url, headers={"User-Agent": "AgriClaw/0.1"})
        with urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())
    except URLError as e:
        print(f"Error fetching weather data: {e}", file=sys.stderr)
        sys.exit(1)

    if "daily" not in data:
        print(f"Unexpected API response: {json.dumps(data, indent=2)}", file=sys.stderr)
        sys.exit(1)

    daily = data["daily"]
    rows = []
    for i, date_str in enumerate(daily["time"]):
        tmax = daily["temperature_2m_max"][i]
        tmin = daily["temperature_2m_min"][i]
        precip = daily["precipitation_sum"][i]
        # Handle missing values
        if tmax is None or tmin is None:
            rows.append({"date": date_str, "tmax": None, "tmin": None, "precip": precip})
        else:
            rows.append({"date": date_str, "tmax": tmax, "tmin": tmin, "precip": precip})

    return rows


def load_weather_csv(filepath):
    """Load user-provided weather CSV. Expected columns: date, tmax, tmin (optional: precip)."""
    rows = []
    with open(filepath, newline="") as f:
        reader = csv.DictReader(f)
        # Normalise column names to lowercase
        for raw_row in reader:
            row = {k.strip().lower(): v.strip() for k, v in raw_row.items()}
            tmax_val = row.get("tmax")
            tmin_val = row.get("tmin")
            precip_val = row.get("precip") or row.get("precipitation") or "0"
            # Auto-detect Fahrenheit: if any value > 60, likely Fahrenheit
            tmax_f = float(tmax_val) if tmax_val else None
            tmin_f = float(tmin_val) if tmin_val else None
            rows.append({
                "date": row["date"],
                "tmax": tmax_f,
                "tmin": tmin_f,
                "precip": float(precip_val) if precip_val else 0.0,
            })

    # Auto-detect Fahrenheit
    valid_temps = [r["tmax"] for r in rows if r["tmax"] is not None]
    if valid_temps and max(valid_temps) > 60:
        print("Detected Fahrenheit — converting to Celsius.", file=sys.stderr)
        for r in rows:
            if r["tmax"] is not None:
                r["tmax"] = (r["tmax"] - 32) * 5 / 9
            if r["tmin"] is not None:
                r["tmin"] = (r["tmin"] - 32) * 5 / 9

    return rows


# ---------------------------------------------------------------------------
# GDD computation
# ---------------------------------------------------------------------------

def compute_gdd(weather_rows, crop):
    """Compute daily and cumulative GDD using modified method with ceiling cap."""
    params = CROP_PARAMS[crop]
    base = params["base"]
    ceiling = params["ceiling"]

    results = []
    cumulative = 0.0
    missing_streak = 0

    for row in weather_rows:
        tmax = row["tmax"]
        tmin = row["tmin"]

        if tmax is None or tmin is None:
            missing_streak += 1
            if missing_streak > 3:
                print(f"Warning: >3 consecutive missing days at {row['date']}. "
                      f"Results after this point may be unreliable.", file=sys.stderr)
            # Interpolate: carry forward previous day's GDD
            daily_gdd = results[-1]["daily_gdd"] if results else 0.0
        else:
            missing_streak = 0
            # Cap Tmax at ceiling
            tmax_capped = min(tmax, ceiling)
            # Floor Tmin at base (no negative contribution)
            tmin_adj = max(tmin, base)
            # Daily GDD
            daily_gdd = max(0.0, (tmax_capped + tmin_adj) / 2.0 - base)

        cumulative += daily_gdd
        frost = tmin is not None and tmin <= 0.0

        results.append({
            "date": row["date"],
            "tmax": tmax,
            "tmin": tmin,
            "precip": row.get("precip", 0.0),
            "daily_gdd": round(daily_gdd, 2),
            "cumulative_gdd": round(cumulative, 2),
            "frost": frost,
        })

    return results


def estimate_stages(gdd_results, crop):
    """Map cumulative GDD to crop growth stages."""
    params = CROP_PARAMS[crop]
    stages = params["stages"]
    reached = []

    for gdd_thresh, stage_code, stage_name in stages:
        # Find first date where cumulative GDD >= threshold
        for row in gdd_results:
            if row["cumulative_gdd"] >= gdd_thresh:
                reached.append({
                    "stage": stage_code,
                    "description": stage_name,
                    "gdd_threshold": gdd_thresh,
                    "date_reached": row["date"],
                    "actual_gdd": row["cumulative_gdd"],
                })
                break
        else:
            reached.append({
                "stage": stage_code,
                "description": stage_name,
                "gdd_threshold": gdd_thresh,
                "date_reached": "Not reached",
                "actual_gdd": gdd_results[-1]["cumulative_gdd"] if gdd_results else 0,
            })

    return reached


def find_frost_dates(gdd_results):
    """Find last spring frost and first fall frost."""
    frost_days = [r for r in gdd_results if r["frost"]]
    if not frost_days:
        return None, None

    # Parse dates to find spring vs fall
    dates_parsed = []
    for r in frost_days:
        dt = datetime.strptime(r["date"], "%Y-%m-%d")
        dates_parsed.append((dt, r["date"]))

    # Last spring frost: last frost before July 1
    year = dates_parsed[0][0].year
    july1 = datetime(year, 7, 1)
    spring_frosts = [(dt, d) for dt, d in dates_parsed if dt < july1]
    last_spring = spring_frosts[-1][1] if spring_frosts else None

    # First fall frost: first frost after July 1
    fall_frosts = [(dt, d) for dt, d in dates_parsed if dt >= july1]
    first_fall = fall_frosts[0][1] if fall_frosts else None

    return last_spring, first_fall


# ---------------------------------------------------------------------------
# Output generation
# ---------------------------------------------------------------------------

def write_csv(filepath, rows, fieldnames):
    """Write a list of dicts to CSV."""
    with open(filepath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(filepath):
    """Compute SHA-256 of a file."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def generate_plots(gdd_results, output_dir):
    """Generate GDD accumulation and daily temperature plots."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
    except ImportError:
        print("matplotlib not installed — skipping plots.", file=sys.stderr)
        return False

    dates = [datetime.strptime(r["date"], "%Y-%m-%d") for r in gdd_results]
    cumulative = [r["cumulative_gdd"] for r in gdd_results]
    tmax = [r["tmax"] for r in gdd_results]
    tmin = [r["tmin"] for r in gdd_results]
    daily_gdd = [r["daily_gdd"] for r in gdd_results]

    fig_dir = os.path.join(output_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    # GDD accumulation plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    ax1.fill_between(dates, cumulative, alpha=0.3, color="#2e7d32")
    ax1.plot(dates, cumulative, color="#2e7d32", linewidth=2)
    ax1.set_ylabel("Cumulative GDD (C-days)")
    ax1.set_title("Growing Degree Day Accumulation")
    ax1.grid(True, alpha=0.3)

    ax2.bar(dates, daily_gdd, color="#66bb6a", alpha=0.7, width=1.0)
    ax2.set_ylabel("Daily GDD (C-days)")
    ax2.set_xlabel("Date")
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax2.xaxis.set_major_locator(mdates.MonthLocator())
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha="right")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "gdd_accumulation.png"), dpi=150)
    plt.close()

    # Daily temperature plot
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.fill_between(dates, tmin, tmax, alpha=0.3, color="#1565c0")
    ax.plot(dates, tmax, color="#e53935", linewidth=1, label="Tmax")
    ax.plot(dates, tmin, color="#1565c0", linewidth=1, label="Tmin")
    ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
    ax.set_ylabel("Temperature (C)")
    ax.set_xlabel("Date")
    ax.set_title("Daily Temperature Range")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "daily_temperature.png"), dpi=150)
    plt.close()

    return True


def generate_report(gdd_results, stages, frost_info, crop, lat, lon, plant_date, output_dir,
                    base_override=None, ceiling_override=None):
    """Generate markdown report."""
    last_spring_frost, first_fall_frost = frost_info
    total_gdd = gdd_results[-1]["cumulative_gdd"]
    total_precip = sum(r["precip"] or 0 for r in gdd_results)
    avg_tmax = sum(r["tmax"] for r in gdd_results if r["tmax"] is not None) / max(1, sum(1 for r in gdd_results if r["tmax"] is not None))
    avg_tmin = sum(r["tmin"] for r in gdd_results if r["tmin"] is not None) / max(1, sum(1 for r in gdd_results if r["tmin"] is not None))
    frost_days = sum(1 for r in gdd_results if r["frost"])
    end_date = gdd_results[-1]["date"]
    params = CROP_PARAMS[crop]

    base_note = " (user override)" if base_override is not None else ""
    ceil_note = " (user override)" if ceiling_override is not None else ""

    lines = [
        f"# GDD Report: {crop.title()}",
        "",
        f"**Location**: {lat}N, {lon}E  ",
        f"**Planting date**: {plant_date}  ",
        f"**Analysis period**: {plant_date} to {end_date}  ",
        f"**Base temperature**: {params['base']}C{base_note}  ",
        f"**Ceiling temperature**: {params['ceiling']}C{ceil_note}  ",
        "",
        "## Season Summary",
        "",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| Total GDD | {total_gdd:.1f} C-days |",
        f"| Total precipitation | {total_precip:.1f} mm |",
        f"| Mean Tmax | {avg_tmax:.1f} C |",
        f"| Mean Tmin | {avg_tmin:.1f} C |",
        f"| Frost days (Tmin <= 0C) | {frost_days} |",
        f"| Last spring frost | {last_spring_frost or 'None'} |",
        f"| First fall frost | {first_fall_frost or 'None'} |",
        "",
        "## Growth Stage Estimates",
        "",
        f"| Stage | Description | GDD Threshold | Date Reached |",
        f"|-------|-------------|---------------|-------------|",
    ]

    for s in stages:
        lines.append(f"| {s['stage']} | {s['description']} | {s['gdd_threshold']} | {s['date_reached']} |")

    # References by crop
    refs = {
        "soybean": [
            "Setiyono et al. (2007). Understanding and modeling the effect of temperature and daylength on soybean phenology under high-yield conditions. Field Crops Research, 100(2-3), 257-271.",
            "Gilmore & Rogers (1958). Heat units as a method of measuring maturity in corn. Agronomy Journal, 50(10), 611-615.",
        ],
        "corn": [
            "Cross & Zuber (1972). Prediction of flowering dates in maize based on different methods of estimating thermal units. Agronomy Journal, 64(3), 351-355.",
            "McMaster & Wilhelm (1997). Growing degree-days: one equation, two interpretations. Agricultural and Forest Meteorology, 87(4), 291-300.",
        ],
        "wheat": [
            "Porter & Gawith (1999). Temperatures and the growth and development of wheat: a review. European Journal of Agronomy, 10(1), 23-36.",
            "McMaster & Smika (1988). Estimation and evaluation of winter wheat phenology in the central Great Plains. Agricultural and Forest Meteorology, 43, 1-18.",
        ],
        "rice": [
            "Yoshida (1981). Fundamentals of Rice Crop Science. IRRI, Los Banos.",
            "Gao et al. (1992). An analysis of the base temperature for computing rice growing degree days. Acta Agronomica Sinica, 18(1), 11-18.",
        ],
    }

    lines += [
        "",
        "## Figures",
        "",
        "![GDD Accumulation](figures/gdd_accumulation.png)",
        "",
        "![Daily Temperature](figures/daily_temperature.png)",
        "",
        "## References",
        "",
        f"Default base ({params['base']}C) and ceiling ({params['ceiling']}C) temperatures for {crop} are based on:",
        "",
    ]
    for ref in refs.get(crop, []):
        lines.append(f"- {ref}")
    lines += [
        "",
        "Growth stage GDD thresholds are approximate mid-maturity values for US Midwest conditions. "
        "Override with `--base` and `--ceiling` for local calibration.",
        "",
        "---",
        f"*Generated by [AgriClaw claw-gdd](https://github.com/bgtamang/AgriClaw) v0.1.0*",
    ]

    report_path = os.path.join(output_dir, "report.md")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def generate_repro_bundle(output_dir, args_dict):
    """Generate reproducibility bundle: commands.sh, environment.yml, checksums.sha256."""
    repro_dir = os.path.join(output_dir, "reproducibility")
    os.makedirs(repro_dir, exist_ok=True)

    # commands.sh
    base_ceil_flags = ""
    if args_dict.get("base") is not None:
        base_ceil_flags += f" --base {args_dict['base']}"
    if args_dict.get("ceiling") is not None:
        base_ceil_flags += f" --ceiling {args_dict['ceiling']}"

    if args_dict.get("demo"):
        cmd = f"python gdd.py --demo{base_ceil_flags}"
    elif args_dict.get("input"):
        cmd = (f"python gdd.py --input {args_dict['input']} "
               f"--crop {args_dict['crop']} --plant-date {args_dict['plant_date']}{base_ceil_flags}")
    else:
        cmd = (f"python gdd.py --lat {args_dict['lat']} --lon {args_dict['lon']} "
               f"--plant-date {args_dict['plant_date']} --crop {args_dict['crop']}{base_ceil_flags}")

    with open(os.path.join(repro_dir, "commands.sh"), "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"# Reproduce this analysis\n")
        f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"pip install pandas matplotlib requests\n")
        f.write(f"{cmd}\n")

    # environment.yml
    with open(os.path.join(repro_dir, "environment.yml"), "w") as f:
        f.write("name: claw-gdd\n")
        f.write("dependencies:\n")
        f.write("  - python>=3.9\n")
        f.write("  - pip:\n")
        f.write("    - pandas\n")
        f.write("    - matplotlib\n")
        f.write("    - requests\n")

    # checksums.sha256
    checksums = []
    for root, dirs, files in os.walk(output_dir):
        for fname in sorted(files):
            fpath = os.path.join(root, fname)
            rel = os.path.relpath(fpath, output_dir)
            if "reproducibility" in rel:
                continue
            checksums.append(f"{sha256_file(fpath)}  {rel}")

    with open(os.path.join(repro_dir, "checksums.sha256"), "w") as f:
        f.write("\n".join(checksums) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="claw-gdd: Growing Degree Day Calculator (AgriClaw)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  python gdd.py --demo\n"
               "  python gdd.py --lat 40.12 --lon -88.24 --plant-date 2025-05-15 --crop soybean\n"
               "  python gdd.py --input weather.csv --crop corn --plant-date 2025-04-20\n",
    )
    parser.add_argument("--demo", action="store_true",
                        help="Run demo: soybean, Champaign IL, 2025 season")
    parser.add_argument("--lat", type=float, help="Latitude (decimal degrees)")
    parser.add_argument("--lon", type=float, help="Longitude (decimal degrees)")
    parser.add_argument("--plant-date", help="Planting date (YYYY-MM-DD)")
    parser.add_argument("--end-date", help="End date (YYYY-MM-DD). Default: 150 days after planting")
    parser.add_argument("--crop", choices=list(CROP_PARAMS.keys()), default="soybean",
                        help="Crop type (default: soybean)")
    parser.add_argument("--base", type=float, default=None,
                        help="Override base temperature (C). Default: crop-specific")
    parser.add_argument("--ceiling", type=float, default=None,
                        help="Override ceiling temperature (C). Default: crop-specific")
    parser.add_argument("--input", help="Path to user weather CSV (date, tmax, tmin columns)")
    parser.add_argument("--output", default="output",
                        help="Output directory (default: output)")
    args = parser.parse_args()

    # Demo mode
    if args.demo:
        args.lat = 40.12
        args.lon = -88.24
        args.plant_date = "2025-05-15"
        args.crop = "soybean"
        args.output = "demo_output"
        print("Running demo: Soybean, Champaign IL (40.12N, 88.24W), planted 2025-05-15")

    # Validate
    if not args.demo and not args.input and (args.lat is None or args.lon is None):
        parser.error("Provide --lat/--lon, --input, or --demo")
    if not args.demo and not args.plant_date:
        parser.error("--plant-date is required")

    plant_date = args.plant_date
    if args.end_date:
        end_date = args.end_date
    else:
        dt = datetime.strptime(plant_date, "%Y-%m-%d") + timedelta(days=150)
        end_date = dt.strftime("%Y-%m-%d")

    # Apply user overrides for base/ceiling
    if args.base is not None:
        CROP_PARAMS[args.crop]["base"] = args.base
    if args.ceiling is not None:
        CROP_PARAMS[args.crop]["ceiling"] = args.ceiling

    # Get weather data
    base_src = "user" if args.base is not None else "default"
    ceil_src = "user" if args.ceiling is not None else "default"
    print(f"Crop: {args.crop} | Base: {CROP_PARAMS[args.crop]['base']}C ({base_src}) | Ceiling: {CROP_PARAMS[args.crop]['ceiling']}C ({ceil_src})")

    if args.input:
        print(f"Loading weather from: {args.input}")
        weather = load_weather_csv(args.input)
    else:
        print(f"Fetching weather for {args.lat}N, {args.lon}E ({plant_date} to {end_date})...")
        weather = fetch_weather(args.lat, args.lon, plant_date, end_date)

    print(f"Got {len(weather)} days of weather data.")

    # Compute GDD
    gdd_results = compute_gdd(weather, args.crop)
    stages = estimate_stages(gdd_results, args.crop)
    frost_info = find_frost_dates(gdd_results)

    # Create output
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    tables_dir = os.path.join(output_dir, "tables")
    os.makedirs(tables_dir, exist_ok=True)

    # Write CSVs
    write_csv(
        os.path.join(tables_dir, "daily_gdd.csv"),
        gdd_results,
        ["date", "tmax", "tmin", "precip", "daily_gdd", "cumulative_gdd", "frost"],
    )
    write_csv(
        os.path.join(tables_dir, "growth_stages.csv"),
        stages,
        ["stage", "description", "gdd_threshold", "date_reached", "actual_gdd"],
    )
    print(f"Wrote daily_gdd.csv ({len(gdd_results)} rows)")

    # Plots
    plots_ok = generate_plots(gdd_results, output_dir)
    if plots_ok:
        print("Generated figures: gdd_accumulation.png, daily_temperature.png")

    # Report
    generate_report(gdd_results, stages, frost_info, args.crop,
                    args.lat or "CSV", args.lon or "input", plant_date, output_dir,
                    base_override=args.base, ceiling_override=args.ceiling)
    print("Generated report.md")

    # Reproducibility bundle
    args_dict = {
        "demo": args.demo, "lat": args.lat, "lon": args.lon,
        "plant_date": plant_date, "crop": args.crop, "input": args.input,
        "base": args.base, "ceiling": args.ceiling,
    }
    generate_repro_bundle(output_dir, args_dict)
    print("Generated reproducibility bundle (commands.sh, environment.yml, checksums.sha256)")

    # Summary
    total_gdd = gdd_results[-1]["cumulative_gdd"]
    print(f"\n{'='*50}")
    print(f"Total GDD accumulated: {total_gdd:.1f} C-days")
    print(f"Growth stages reached:")
    for s in stages:
        status = s["date_reached"] if s["date_reached"] != "Not reached" else "— not reached"
        print(f"  {s['stage']:12s} ({s['gdd_threshold']:>5d} GDD) -> {status}")
    last_spring, first_fall = frost_info
    if last_spring:
        print(f"Last spring frost: {last_spring}")
    if first_fall:
        print(f"First fall frost: {first_fall}")
    print(f"\nOutput: {os.path.abspath(output_dir)}/")


if __name__ == "__main__":
    main()
