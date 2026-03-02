#!/usr/bin/env python3
"""
claw-weather: Weather Season Summary for AgriClaw

Generates a comprehensive weather season summary for a field location using
the Open-Meteo archive API (free, no key needed). Includes temperature and
precipitation analysis, drought indices, dry spell detection, extreme event
counts, and a full reproducibility bundle.

Usage:
    python weather.py --demo
    python weather.py --lat 40.12 --lon -88.24 --start-date 2025-04-01 --end-date 2025-10-31
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
# Weather data fetching (Open-Meteo archive, free, no API key)
# ---------------------------------------------------------------------------

def fetch_weather(lat, lon, start_date, end_date):
    """Fetch daily weather variables from Open-Meteo historical archive."""
    url = (
        f"https://archive-api.open-meteo.com/v1/archive?"
        f"latitude={lat}&longitude={lon}"
        f"&start_date={start_date}&end_date={end_date}"
        f"&daily=temperature_2m_max,temperature_2m_min,precipitation_sum,"
        f"et0_fao_evapotranspiration,rain_sum,shortwave_radiation_sum"
        f"&timezone=auto"
    )
    try:
        req = Request(url, headers={"User-Agent": "AgriClaw/0.1"})
        with urlopen(req, timeout=60) as resp:
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
        rows.append({
            "date": date_str,
            "tmax": daily["temperature_2m_max"][i],
            "tmin": daily["temperature_2m_min"][i],
            "precip": daily["precipitation_sum"][i],
            "et0": daily["et0_fao_evapotranspiration"][i],
            "rain": daily["rain_sum"][i],
            "radiation": daily["shortwave_radiation_sum"][i],
        })

    return rows


# ---------------------------------------------------------------------------
# Analysis functions
# ---------------------------------------------------------------------------

def safe_val(v, default=0.0):
    """Return v if not None, else default."""
    return v if v is not None else default


def compute_season_summary(weather_rows):
    """Compute season-level summary statistics."""
    n = len(weather_rows)
    valid_tmax = [r["tmax"] for r in weather_rows if r["tmax"] is not None]
    valid_tmin = [r["tmin"] for r in weather_rows if r["tmin"] is not None]

    total_precip = sum(safe_val(r["precip"]) for r in weather_rows)
    total_et0 = sum(safe_val(r["et0"]) for r in weather_rows)
    total_radiation = sum(safe_val(r["radiation"]) for r in weather_rows)

    mean_tmax = sum(valid_tmax) / len(valid_tmax) if valid_tmax else 0.0
    mean_tmin = sum(valid_tmin) / len(valid_tmin) if valid_tmin else 0.0
    max_tmax = max(valid_tmax) if valid_tmax else 0.0
    min_tmin = min(valid_tmin) if valid_tmin else 0.0

    return {
        "n_days": n,
        "total_precip_mm": round(total_precip, 1),
        "total_et0_mm": round(total_et0, 1),
        "total_radiation_mj_m2": round(total_radiation, 1),
        "mean_tmax_c": round(mean_tmax, 1),
        "mean_tmin_c": round(mean_tmin, 1),
        "max_tmax_c": round(max_tmax, 1),
        "min_tmin_c": round(min_tmin, 1),
        "mean_temp_c": round((mean_tmax + mean_tmin) / 2, 1),
    }


def compute_monthly_breakdown(weather_rows):
    """Compute monthly summary table."""
    months = {}
    for r in weather_rows:
        ym = r["date"][:7]  # YYYY-MM
        if ym not in months:
            months[ym] = []
        months[ym].append(r)

    results = []
    for ym in sorted(months.keys()):
        rows = months[ym]
        n = len(rows)
        valid_tmax = [r["tmax"] for r in rows if r["tmax"] is not None]
        valid_tmin = [r["tmin"] for r in rows if r["tmin"] is not None]

        total_precip = sum(safe_val(r["precip"]) for r in rows)
        total_et0 = sum(safe_val(r["et0"]) for r in rows)
        total_radiation = sum(safe_val(r["radiation"]) for r in rows)
        mean_tmax = sum(valid_tmax) / len(valid_tmax) if valid_tmax else 0.0
        mean_tmin = sum(valid_tmin) / len(valid_tmin) if valid_tmin else 0.0
        deficit = total_precip - total_et0

        results.append({
            "month": ym,
            "n_days": n,
            "mean_tmax_c": round(mean_tmax, 1),
            "mean_tmin_c": round(mean_tmin, 1),
            "total_precip_mm": round(total_precip, 1),
            "total_et0_mm": round(total_et0, 1),
            "precip_deficit_mm": round(deficit, 1),
            "total_radiation_mj_m2": round(total_radiation, 1),
        })

    return results


def compute_drought_index(monthly_data):
    """Compute simple precipitation deficit (P - ET0) by month.

    Negative values indicate potential water stress. This is a simplified
    climatic water balance following the approach described in Allen et al.
    (1998) FAO Irrigation and Drainage Paper 56, Section 2.
    """
    # Already computed in monthly breakdown as precip_deficit_mm
    return monthly_data


def detect_dry_spells(weather_rows, threshold_mm=1.0):
    """Detect consecutive dry days (precip < threshold).

    Returns list of dry spells sorted by length (longest first).
    """
    spells = []
    current_start = None
    current_length = 0

    for r in weather_rows:
        precip = safe_val(r["precip"])
        if precip < threshold_mm:
            if current_start is None:
                current_start = r["date"]
            current_length += 1
        else:
            if current_length > 0:
                spells.append({
                    "start": current_start,
                    "end": weather_rows[weather_rows.index(r) - 1]["date"] if current_length > 0 else current_start,
                    "length_days": current_length,
                })
            current_start = None
            current_length = 0

    # Handle spell running to end of data
    if current_length > 0:
        spells.append({
            "start": current_start,
            "end": weather_rows[-1]["date"],
            "length_days": current_length,
        })

    spells.sort(key=lambda x: x["length_days"], reverse=True)
    return spells


def count_extreme_events(weather_rows):
    """Count extreme weather events.

    - Days with Tmax > 35C (heat stress)
    - Days with Tmin < 0C (frost)
    - Days with precip > 25mm (heavy rain)
    """
    hot_days = sum(1 for r in weather_rows if r["tmax"] is not None and r["tmax"] > 35.0)
    frost_days = sum(1 for r in weather_rows if r["tmin"] is not None and r["tmin"] < 0.0)
    heavy_rain_days = sum(1 for r in weather_rows if r["precip"] is not None and r["precip"] > 25.0)

    hot_dates = [r["date"] for r in weather_rows if r["tmax"] is not None and r["tmax"] > 35.0]
    frost_dates = [r["date"] for r in weather_rows if r["tmin"] is not None and r["tmin"] < 0.0]
    heavy_rain_dates = [r["date"] for r in weather_rows if r["precip"] is not None and r["precip"] > 25.0]

    return {
        "hot_days_above_35c": hot_days,
        "frost_days_below_0c": frost_days,
        "heavy_rain_days_above_25mm": heavy_rain_days,
        "hot_dates": hot_dates,
        "frost_dates": frost_dates,
        "heavy_rain_dates": heavy_rain_dates,
    }


# ---------------------------------------------------------------------------
# Output generation
# ---------------------------------------------------------------------------

def write_csv(filepath, rows, fieldnames):
    """Write a list of dicts to CSV."""
    with open(filepath, "w", newline="", encoding="utf-8") as f:
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


def generate_plots(weather_rows, monthly_data, output_dir):
    """Generate weather summary figures."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
    except ImportError:
        print("matplotlib not installed -- skipping plots.", file=sys.stderr)
        return False

    dates = [datetime.strptime(r["date"], "%Y-%m-%d") for r in weather_rows]
    tmax = [r["tmax"] for r in weather_rows]
    tmin = [r["tmin"] for r in weather_rows]
    precip = [safe_val(r["precip"]) for r in weather_rows]
    et0 = [safe_val(r["et0"]) for r in weather_rows]
    radiation = [safe_val(r["radiation"]) for r in weather_rows]

    fig_dir = os.path.join(output_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    # -----------------------------------------------------------------------
    # Figure 1: Temperature + Precipitation combo chart
    # -----------------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(12, 5))

    # Precipitation as bars on secondary axis
    ax2 = ax1.twinx()
    ax2.bar(dates, precip, color="#4fc3f7", alpha=0.5, width=1.0, label="Precip (mm)")
    ax2.set_ylabel("Precipitation (mm)")
    ax2.set_ylim(0, max(precip) * 2.5 if max(precip) > 0 else 10)

    # Temperature as lines on primary axis
    ax1.fill_between(dates, tmin, tmax, alpha=0.2, color="#ef5350")
    ax1.plot(dates, tmax, color="#e53935", linewidth=1.2, label="Tmax")
    ax1.plot(dates, tmin, color="#1565c0", linewidth=1.2, label="Tmin")
    ax1.set_ylabel("Temperature (C)")
    ax1.set_xlabel("Date")
    ax1.set_title("Daily Temperature and Precipitation")
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax1.xaxis.set_major_locator(mdates.MonthLocator())
    plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha="right")
    ax1.grid(True, alpha=0.3)

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left")

    ax1.set_zorder(ax2.get_zorder() + 1)
    ax1.patch.set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "temp_precip.png"), dpi=150)
    plt.close()

    # -----------------------------------------------------------------------
    # Figure 2: Cumulative precipitation vs cumulative ET0 (water balance)
    # -----------------------------------------------------------------------
    cum_precip = []
    cum_et0 = []
    running_p = 0.0
    running_e = 0.0
    for r in weather_rows:
        running_p += safe_val(r["precip"])
        running_e += safe_val(r["et0"])
        cum_precip.append(running_p)
        cum_et0.append(running_e)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(dates, cum_precip, color="#1565c0", linewidth=2, label="Cumulative P")
    ax.plot(dates, cum_et0, color="#e53935", linewidth=2, label="Cumulative ET0")
    ax.fill_between(dates, cum_precip, cum_et0, where=[p >= e for p, e in zip(cum_precip, cum_et0)],
                    alpha=0.2, color="#4caf50", label="Surplus (P > ET0)")
    ax.fill_between(dates, cum_precip, cum_et0, where=[p < e for p, e in zip(cum_precip, cum_et0)],
                    alpha=0.2, color="#ef5350", label="Deficit (P < ET0)")
    ax.set_ylabel("Cumulative (mm)")
    ax.set_xlabel("Date")
    ax.set_title("Cumulative Precipitation vs. Reference ET0 (Water Balance)")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "water_balance.png"), dpi=150)
    plt.close()

    # -----------------------------------------------------------------------
    # Figure 3: Monthly summary bar chart
    # -----------------------------------------------------------------------
    month_labels = [m["month"] for m in monthly_data]
    month_precip = [m["total_precip_mm"] for m in monthly_data]
    month_et0 = [m["total_et0_mm"] for m in monthly_data]
    month_tmax = [m["mean_tmax_c"] for m in monthly_data]
    month_tmin = [m["mean_tmin_c"] for m in monthly_data]

    x_pos = list(range(len(month_labels)))
    bar_width = 0.35

    fig, ax1 = plt.subplots(figsize=(10, 5))

    bars1 = ax1.bar([x - bar_width / 2 for x in x_pos], month_precip,
                    bar_width, label="Precip (mm)", color="#4fc3f7")
    bars2 = ax1.bar([x + bar_width / 2 for x in x_pos], month_et0,
                    bar_width, label="ET0 (mm)", color="#ef9a9a")
    ax1.set_ylabel("Total (mm)")
    ax1.set_xlabel("Month")
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(month_labels, rotation=45, ha="right")

    # Temperature on secondary axis
    ax2 = ax1.twinx()
    ax2.plot(x_pos, month_tmax, color="#e53935", marker="o", linewidth=2, label="Mean Tmax")
    ax2.plot(x_pos, month_tmin, color="#1565c0", marker="s", linewidth=2, label="Mean Tmin")
    ax2.set_ylabel("Temperature (C)")

    ax1.set_title("Monthly Weather Summary")

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left")
    ax1.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "monthly_summary.png"), dpi=150)
    plt.close()

    return True


def generate_report(weather_rows, summary, monthly_data, dry_spells, extremes,
                    lat, lon, start_date, end_date, output_dir):
    """Generate markdown report."""

    lines = [
        "# Weather Season Summary",
        "",
        f"**Location**: {lat}\u00b0N, {lon}\u00b0E  ",
        f"**Period**: {start_date} to {end_date}  ",
        f"**Days analysed**: {summary['n_days']}  ",
        "",
        "## Season Summary",
        "",
        "| Metric | Value |",
        "|--------|-------|",
        f"| Total precipitation | {summary['total_precip_mm']} mm |",
        f"| Total reference ET0 | {summary['total_et0_mm']} mm |",
        f"| Season water balance (P - ET0) | {round(summary['total_precip_mm'] - summary['total_et0_mm'], 1)} mm |",
        f"| Total shortwave radiation | {summary['total_radiation_mj_m2']} MJ/m\u00b2 |",
        f"| Mean Tmax | {summary['mean_tmax_c']}\u00b0C |",
        f"| Mean Tmin | {summary['mean_tmin_c']}\u00b0C |",
        f"| Mean temperature | {summary['mean_temp_c']}\u00b0C |",
        f"| Season max temperature | {summary['max_tmax_c']}\u00b0C |",
        f"| Season min temperature | {summary['min_tmin_c']}\u00b0C |",
        "",
        "## Monthly Breakdown",
        "",
        "| Month | Days | Mean Tmax (\u00b0C) | Mean Tmin (\u00b0C) | Precip (mm) | ET0 (mm) | P - ET0 (mm) | Radiation (MJ/m\u00b2) |",
        "|-------|------|--------------|--------------|-------------|----------|---------------|---------------------|",
    ]

    for m in monthly_data:
        lines.append(
            f"| {m['month']} | {m['n_days']} | {m['mean_tmax_c']} | {m['mean_tmin_c']} "
            f"| {m['total_precip_mm']} | {m['total_et0_mm']} | {m['precip_deficit_mm']} "
            f"| {m['total_radiation_mj_m2']} |"
        )

    # Drought index section
    lines += [
        "",
        "## Drought Index (Precipitation Deficit)",
        "",
        "Simple climatic water balance: P - ET0 by month. Negative values indicate",
        "months where atmospheric evaporative demand exceeded precipitation, suggesting",
        "potential crop water stress without irrigation.",
        "",
        "| Month | P (mm) | ET0 (mm) | P - ET0 (mm) | Status |",
        "|-------|--------|----------|---------------|--------|",
    ]

    for m in monthly_data:
        deficit = m["precip_deficit_mm"]
        if deficit >= 0:
            status = "Surplus"
        elif deficit > -50:
            status = "Mild deficit"
        elif deficit > -100:
            status = "Moderate deficit"
        else:
            status = "Severe deficit"
        lines.append(
            f"| {m['month']} | {m['total_precip_mm']} | {m['total_et0_mm']} "
            f"| {deficit} | {status} |"
        )

    # Dry spells
    lines += [
        "",
        "## Dry Spell Analysis",
        "",
        "Consecutive days with precipitation < 1 mm.",
        "",
    ]
    if dry_spells:
        longest = dry_spells[0]
        lines.append(
            f"**Longest dry spell**: {longest['length_days']} days "
            f"({longest['start']} to {longest['end']})"
        )
        lines += [
            "",
            "| Rank | Start | End | Length (days) |",
            "|------|-------|-----|---------------|",
        ]
        for i, spell in enumerate(dry_spells[:10]):
            lines.append(
                f"| {i + 1} | {spell['start']} | {spell['end']} | {spell['length_days']} |"
            )
    else:
        lines.append("No dry spells detected (every day had >= 1 mm precipitation).")

    # Extreme events
    lines += [
        "",
        "## Extreme Events",
        "",
        "| Event | Count |",
        "|-------|-------|",
        f"| Days with Tmax > 35\u00b0C (heat stress) | {extremes['hot_days_above_35c']} |",
        f"| Days with Tmin < 0\u00b0C (frost) | {extremes['frost_days_below_0c']} |",
        f"| Days with precip > 25 mm (heavy rain) | {extremes['heavy_rain_days_above_25mm']} |",
        "",
    ]

    if extremes["hot_dates"]:
        lines.append(f"**Heat stress dates**: {', '.join(extremes['hot_dates'])}")
        lines.append("")
    if extremes["frost_dates"]:
        lines.append(f"**Frost dates**: {', '.join(extremes['frost_dates'])}")
        lines.append("")
    if extremes["heavy_rain_dates"]:
        lines.append(f"**Heavy rain dates**: {', '.join(extremes['heavy_rain_dates'])}")
        lines.append("")

    # Figures
    lines += [
        "## Figures",
        "",
        "![Temperature and Precipitation](figures/temp_precip.png)",
        "",
        "![Water Balance](figures/water_balance.png)",
        "",
        "![Monthly Summary](figures/monthly_summary.png)",
        "",
    ]

    # References
    lines += [
        "## References",
        "",
        "- Allen, R.G., Pereira, L.S., Raes, D., & Smith, M. (1998). Crop evapotranspiration: "
        "Guidelines for computing crop water requirements. FAO Irrigation and Drainage Paper 56. Rome: FAO.",
        "- Doorenbos, J. & Pruitt, W.O. (1977). Guidelines for predicting crop water requirements. "
        "FAO Irrigation and Drainage Paper 24. Rome: FAO.",
        "- Open-Meteo (2023). Free Weather API. https://open-meteo.com/",
        "",
        "ET0 values from Open-Meteo use the FAO Penman-Monteith equation (Allen et al., 1998).",
        "The precipitation deficit (P - ET0) is a simplified climatic water balance indicator;",
        "it does not account for soil water storage, runoff, or deep percolation.",
        "",
        "---",
        f"*Generated by [AgriClaw claw-weather](https://github.com/bgtamang/AgriClaw) v0.1.0*",
    ]

    report_path = os.path.join(output_dir, "report.md")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def generate_repro_bundle(output_dir, args_dict):
    """Generate reproducibility bundle: commands.sh, environment.yml, checksums.sha256."""
    repro_dir = os.path.join(output_dir, "reproducibility")
    os.makedirs(repro_dir, exist_ok=True)

    # commands.sh
    if args_dict.get("demo"):
        cmd = "python weather.py --demo"
    else:
        cmd = (
            f"python weather.py --lat {args_dict['lat']} --lon {args_dict['lon']} "
            f"--start-date {args_dict['start_date']} --end-date {args_dict['end_date']}"
        )

    with open(os.path.join(repro_dir, "commands.sh"), "w", encoding="utf-8") as f:
        f.write("#!/bin/bash\n")
        f.write(f"# Reproduce this analysis\n")
        f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"pip install matplotlib\n")
        f.write(f"{cmd}\n")

    # environment.yml
    with open(os.path.join(repro_dir, "environment.yml"), "w", encoding="utf-8") as f:
        f.write("name: claw-weather\n")
        f.write("dependencies:\n")
        f.write("  - python>=3.9\n")
        f.write("  - pip:\n")
        f.write("    - matplotlib\n")

    # checksums.sha256
    checksums = []
    for root, dirs, files in os.walk(output_dir):
        for fname in sorted(files):
            fpath = os.path.join(root, fname)
            rel = os.path.relpath(fpath, output_dir)
            if "reproducibility" in rel:
                continue
            checksums.append(f"{sha256_file(fpath)}  {rel}")

    with open(os.path.join(repro_dir, "checksums.sha256"), "w", encoding="utf-8") as f:
        f.write("\n".join(checksums) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="claw-weather: Weather Season Summary (AgriClaw)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  python weather.py --demo\n"
               "  python weather.py --lat 40.12 --lon -88.24 --start-date 2025-04-01 --end-date 2025-10-31\n",
    )
    parser.add_argument("--demo", action="store_true",
                        help="Run demo: Champaign IL, 2025 growing season (Apr-Oct)")
    parser.add_argument("--lat", type=float, help="Latitude (decimal degrees)")
    parser.add_argument("--lon", type=float, help="Longitude (decimal degrees)")
    parser.add_argument("--start-date", help="Start date (YYYY-MM-DD)")
    parser.add_argument("--end-date", help="End date (YYYY-MM-DD)")
    parser.add_argument("--output", default="output",
                        help="Output directory (default: output)")
    args = parser.parse_args()

    # Demo mode
    if args.demo:
        args.lat = 40.12
        args.lon = -88.24
        args.start_date = "2025-04-01"
        args.end_date = "2025-10-31"
        args.output = "demo_output"
        print("Running demo: Champaign IL (40.12N, 88.24W), 2025 growing season (Apr 1 - Oct 31)")

    # Validate
    if not args.demo:
        if args.lat is None or args.lon is None:
            parser.error("Provide --lat and --lon, or use --demo")
        if args.start_date is None or args.end_date is None:
            parser.error("Provide --start-date and --end-date, or use --demo")

    lat = args.lat
    lon = args.lon
    start_date = args.start_date
    end_date = args.end_date

    print(f"Fetching weather for {lat}N, {lon}E ({start_date} to {end_date})...")
    weather_rows = fetch_weather(lat, lon, start_date, end_date)
    print(f"Got {len(weather_rows)} days of weather data.")

    # Analysis
    print("Computing season summary...")
    summary = compute_season_summary(weather_rows)
    monthly_data = compute_monthly_breakdown(weather_rows)
    dry_spells = detect_dry_spells(weather_rows, threshold_mm=1.0)
    extremes = count_extreme_events(weather_rows)

    # Create output directories
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    tables_dir = os.path.join(output_dir, "tables")
    os.makedirs(tables_dir, exist_ok=True)

    # Write CSVs
    write_csv(
        os.path.join(tables_dir, "daily_weather.csv"),
        weather_rows,
        ["date", "tmax", "tmin", "precip", "et0", "rain", "radiation"],
    )
    print(f"Wrote daily_weather.csv ({len(weather_rows)} rows)")

    # Monthly breakdown CSV (exclude None-causing fields for extremes)
    write_csv(
        os.path.join(tables_dir, "monthly_summary.csv"),
        monthly_data,
        ["month", "n_days", "mean_tmax_c", "mean_tmin_c",
         "total_precip_mm", "total_et0_mm", "precip_deficit_mm", "total_radiation_mj_m2"],
    )
    print(f"Wrote monthly_summary.csv ({len(monthly_data)} rows)")

    # Dry spells CSV
    if dry_spells:
        write_csv(
            os.path.join(tables_dir, "dry_spells.csv"),
            dry_spells,
            ["start", "end", "length_days"],
        )
        print(f"Wrote dry_spells.csv ({len(dry_spells)} spells)")

    # Extreme events CSV
    extreme_rows = []
    for d in extremes["hot_dates"]:
        extreme_rows.append({"date": d, "event": "tmax_above_35c"})
    for d in extremes["frost_dates"]:
        extreme_rows.append({"date": d, "event": "tmin_below_0c"})
    for d in extremes["heavy_rain_dates"]:
        extreme_rows.append({"date": d, "event": "precip_above_25mm"})
    if extreme_rows:
        extreme_rows.sort(key=lambda x: x["date"])
        write_csv(
            os.path.join(tables_dir, "extreme_events.csv"),
            extreme_rows,
            ["date", "event"],
        )
        print(f"Wrote extreme_events.csv ({len(extreme_rows)} events)")

    # Plots
    plots_ok = generate_plots(weather_rows, monthly_data, output_dir)
    if plots_ok:
        print("Generated figures: temp_precip.png, water_balance.png, monthly_summary.png")

    # Report
    generate_report(weather_rows, summary, monthly_data, dry_spells, extremes,
                    lat, lon, start_date, end_date, output_dir)
    print("Generated report.md")

    # Reproducibility bundle
    args_dict = {
        "demo": args.demo,
        "lat": lat,
        "lon": lon,
        "start_date": start_date,
        "end_date": end_date,
    }
    generate_repro_bundle(output_dir, args_dict)
    print("Generated reproducibility bundle (commands.sh, environment.yml, checksums.sha256)")

    # Summary to stdout
    print("")
    print("=" * 60)
    print("WEATHER SEASON SUMMARY")
    print("=" * 60)
    print(f"  Location:           {lat}N, {lon}E")
    print(f"  Period:             {start_date} to {end_date} ({summary['n_days']} days)")
    print(f"  Total precip:       {summary['total_precip_mm']} mm")
    print(f"  Total ET0:          {summary['total_et0_mm']} mm")
    print(f"  Water balance:      {round(summary['total_precip_mm'] - summary['total_et0_mm'], 1)} mm")
    print(f"  Mean Tmax:          {summary['mean_tmax_c']} C")
    print(f"  Mean Tmin:          {summary['mean_tmin_c']} C")
    print(f"  Season max temp:    {summary['max_tmax_c']} C")
    print(f"  Season min temp:    {summary['min_tmin_c']} C")
    print(f"  Total radiation:    {summary['total_radiation_mj_m2']} MJ/m2")
    print(f"  Heat days (>35C):   {extremes['hot_days_above_35c']}")
    print(f"  Frost days (<0C):   {extremes['frost_days_below_0c']}")
    print(f"  Heavy rain (>25mm): {extremes['heavy_rain_days_above_25mm']}")
    if dry_spells:
        longest = dry_spells[0]
        print(f"  Longest dry spell:  {longest['length_days']} days ({longest['start']} to {longest['end']})")
    else:
        print(f"  Longest dry spell:  0 days")
    print(f"\nOutput: {os.path.abspath(output_dir)}/")


if __name__ == "__main__":
    main()
