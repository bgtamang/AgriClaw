#!/usr/bin/env python3
"""
claw-soil: Soil Property Lookup via USDA Soil Data Access (SDA)

Queries the SSURGO database through the SDA REST API to retrieve soil
properties for a given location: texture, organic matter, pH, drainage,
taxonomic classification, available water capacity, and Ksat.

Usage:
    python soil.py --demo
    python soil.py --lat 40.08 --lon -88.20
    python soil.py --lat 42.03 --lon -93.47 --output my_soil_report
"""

import argparse
import csv
import hashlib
import json
import os
import sys
from datetime import datetime
from urllib.request import urlopen, Request
from urllib.error import URLError

# ---------------------------------------------------------------------------
# USDA Soil Data Access (SDA) API
# https://sdmdataaccess.sc.egov.usda.gov/
# ---------------------------------------------------------------------------

SDA_URL = "https://sdmdataaccess.sc.egov.usda.gov/tabular/post.rest"


def sda_query(sql):
    """Execute a SQL query against the USDA Soil Data Access REST API.

    Returns a list of dicts (one per row), or an empty list on failure.
    The SDA API returns JSON with a 'Table' key containing rows.
    """
    payload = json.dumps({"query": sql, "format": "JSON+COLUMNNAME"})
    req = Request(
        SDA_URL,
        data=payload.encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "User-Agent": "AgriClaw/0.1",
        },
    )
    try:
        with urlopen(req, timeout=90) as resp:
            raw = resp.read().decode("utf-8")
    except (URLError, TimeoutError, OSError) as e:
        print("Error contacting USDA Soil Data Access: %s" % e, file=sys.stderr)
        return None

    try:
        data = json.loads(raw)
    except json.JSONDecodeError:
        print("SDA returned non-JSON response. First 500 chars:", file=sys.stderr)
        print(raw[:500], file=sys.stderr)
        return None

    if "Table" not in data:
        # Empty dict {} means no results (not an error); other shapes may be errors
        if data == {}:
            return []
        print("SDA response has no 'Table' key: %s" % json.dumps(data, indent=2)[:500],
              file=sys.stderr)
        return None

    rows = data["Table"]
    if len(rows) < 2:
        # First row is column names, need at least one data row
        return []

    columns = rows[0]
    result = []
    for row in rows[1:]:
        result.append(dict(zip(columns, row)))
    return result


# ---------------------------------------------------------------------------
# Soil data retrieval (two-step approach for reliability)
# ---------------------------------------------------------------------------

def get_mukey_for_point(lat, lon):
    """Get the map unit key (mukey) for a geographic point."""
    sql = (
        "SELECT * FROM SDA_Get_Mukey_from_intersection_with_WktWgs84("
        "'POINT(%s %s)')" % (lon, lat)
    )
    rows = sda_query(sql)
    if rows is None or len(rows) == 0:
        return None
    # The function returns mukey and other fields
    mukey = rows[0].get("mukey") or rows[0].get("MUKEY")
    if mukey is None:
        # Try first value in the dict
        mukey = list(rows[0].values())[0]
    return str(mukey)


def get_soil_properties(mukey):
    """Get soil properties for a map unit key, including all horizons of major components."""
    sql = """
    SELECT
        mapunit.muname,
        component.compname,
        component.comppct_r,
        component.taxclname,
        component.taxorder,
        component.taxsuborder,
        component.drainagecl,
        component.majcompflag,
        chorizon.hzdept_r,
        chorizon.hzdepb_r,
        chorizon.hzname,
        chorizon.sandtotal_r,
        chorizon.silttotal_r,
        chorizon.claytotal_r,
        chorizon.om_r,
        chorizon.ph1to1h2o_r,
        chorizon.awc_r,
        chorizon.ksat_r,
        chorizon.dbthirdbar_r,
        chorizon.cec7_r
    FROM mapunit
    INNER JOIN component ON component.mukey = mapunit.mukey
    INNER JOIN chorizon ON chorizon.cokey = component.cokey
    WHERE mapunit.mukey = '%s'
      AND component.majcompflag = 'Yes'
    ORDER BY component.comppct_r DESC, chorizon.hzdept_r ASC
    """ % mukey
    return sda_query(sql)


def get_mapunit_info(mukey):
    """Get map unit level information."""
    sql = """
    SELECT
        mapunit.muname,
        mapunit.mukind,
        mapunit.muacres,
        legend.areaname,
        legend.areasymbol
    FROM mapunit
    INNER JOIN legend ON legend.lkey = mapunit.lkey
    WHERE mapunit.mukey = '%s'
    """ % mukey
    return sda_query(sql)


# ---------------------------------------------------------------------------
# Soil texture triangle classification
# ---------------------------------------------------------------------------

def classify_texture(sand, silt, clay):
    """Classify soil texture using the USDA soil texture triangle.

    Parameters are percentages (0-100). Returns texture class name.

    Reference:
        USDA-NRCS Soil Survey Manual, Chapter 3.
        Soil Science Division Staff (2017). Soil Survey Manual.
        USDA Handbook 18, Washington, D.C.
    """
    if sand is None or silt is None or clay is None:
        return "Unknown"

    sand = float(sand)
    silt = float(silt)
    clay = float(clay)

    # Validate percentages sum to ~100
    total = sand + silt + clay
    if total < 90 or total > 110:
        return "Unknown (invalid %)"

    # USDA texture classes (checked against official NRCS triangle)
    if clay >= 40 and sand <= 45 and silt < 40:
        return "Clay"
    elif clay >= 40 and silt >= 40:
        return "Silty clay"
    elif clay >= 35 and sand > 45:
        return "Sandy clay"
    elif clay >= 27 and clay < 40 and sand > 20 and sand <= 45:
        return "Clay loam"
    elif clay >= 27 and clay < 40 and sand <= 20:
        return "Silty clay loam"
    elif clay >= 20 and clay < 35 and silt < 28 and sand > 45:
        return "Sandy clay loam"
    elif silt >= 80 and clay < 12:
        return "Silt"
    elif silt >= 50 and clay >= 12 and clay < 27:
        return "Silty clay loam"
    elif silt >= 50 and clay < 12:
        return "Silt loam"
    elif silt >= 28 and silt < 50 and clay < 27 and sand <= 52:
        return "Loam"
    elif clay < 7 and silt < 50 and sand > 52:
        if sand > 85:
            return "Sand"
        else:
            return "Loamy sand"
    elif clay < 20 and silt < 50 and sand >= 52:
        return "Sandy loam"
    elif silt >= 50:
        return "Silt loam"
    else:
        return "Loam"


# ---------------------------------------------------------------------------
# Interpret soil properties for agronomic use
# ---------------------------------------------------------------------------

def interpret_drainage(drainage_class):
    """Provide agronomic interpretation of drainage class."""
    if drainage_class is None:
        return "No data"
    dc = str(drainage_class).lower()
    interpretations = {
        "excessively drained": "Very low water-holding capacity; may need irrigation. Risk of nutrient leaching.",
        "somewhat excessively drained": "Low water retention; consider irrigation for high-demand crops.",
        "well drained": "Ideal for most crops. Good aeration and root development.",
        "moderately well drained": "Generally suitable; may have brief seasonal wetness at depth.",
        "somewhat poorly drained": "Seasonal high water table; may need tile drainage for row crops.",
        "poorly drained": "High water table; tile drainage typically required for row crops.",
        "very poorly drained": "Standing water common; limited to wetland species unless heavily drained.",
    }
    return interpretations.get(dc, "See NRCS drainage class guide.")


def interpret_ph(ph):
    """Provide agronomic interpretation of soil pH."""
    if ph is None:
        return "No data"
    ph = float(ph)
    if ph < 4.5:
        return "Extremely acidic -- severe Al toxicity risk; lime heavily"
    elif ph < 5.5:
        return "Strongly acidic -- lime recommended for most crops"
    elif ph < 6.0:
        return "Moderately acidic -- lime may benefit sensitive crops"
    elif ph < 6.5:
        return "Slightly acidic -- suitable for most crops"
    elif ph < 7.3:
        return "Neutral -- optimal for most crops"
    elif ph < 7.8:
        return "Slightly alkaline -- monitor micronutrient availability"
    elif ph < 8.4:
        return "Moderately alkaline -- Fe/Zn deficiency risk"
    else:
        return "Strongly alkaline -- sodic conditions likely"


def interpret_om(om_pct):
    """Provide agronomic interpretation of organic matter percentage."""
    if om_pct is None:
        return "No data"
    om = float(om_pct)
    if om < 1.0:
        return "Very low -- limited nutrient supply and water retention"
    elif om < 2.0:
        return "Low -- consider cover crops and organic amendments"
    elif om < 4.0:
        return "Moderate -- typical for Midwest cropland"
    elif om < 6.0:
        return "High -- good nutrient supply and water retention"
    else:
        return "Very high -- excellent soil health indicator"


# ---------------------------------------------------------------------------
# Output generation
# ---------------------------------------------------------------------------

def write_csv(filepath, rows, fieldnames):
    """Write a list of dicts to CSV."""
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(filepath):
    """Compute SHA-256 of a file."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def safe_float(val, fmt="%.1f"):
    """Safely format a value as float, returning '--' if None or empty."""
    if val is None or str(val).strip() == "":
        return "--"
    try:
        return fmt % float(val)
    except (ValueError, TypeError):
        return "--"


def generate_profile_plot(horizons, output_dir):
    """Generate a soil profile diagram showing horizons with depth and properties."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    except ImportError:
        print("matplotlib not installed -- skipping plots.", file=sys.stderr)
        return False

    fig_dir = os.path.join(output_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    if not horizons:
        return False

    # Build horizon data for the dominant component
    comp_name = horizons[0].get("compname", "Unknown")
    profile_data = []
    for h in horizons:
        if h.get("compname") != comp_name:
            break  # Only plot dominant component
        top = h.get("hzdept_r")
        bot = h.get("hzdepb_r")
        if top is None or bot is None:
            continue
        top = float(top)
        bot = float(bot)
        sand = float(h["sandtotal_r"]) if h.get("sandtotal_r") else 0
        silt = float(h["silttotal_r"]) if h.get("silttotal_r") else 0
        clay = float(h["claytotal_r"]) if h.get("claytotal_r") else 0
        om = float(h["om_r"]) if h.get("om_r") else 0
        ph = float(h["ph1to1h2o_r"]) if h.get("ph1to1h2o_r") else 0
        texture = classify_texture(
            h.get("sandtotal_r"), h.get("silttotal_r"), h.get("claytotal_r")
        )
        hz_name = h.get("hzname") or ""
        profile_data.append({
            "top": top, "bot": bot,
            "sand": sand, "silt": silt, "clay": clay,
            "om": om, "ph": ph,
            "texture": texture, "hzname": hz_name,
        })

    if not profile_data:
        return False

    max_depth = max(h["bot"] for h in profile_data)

    # Color based on organic matter (darker = more OM)
    def om_color(om_val):
        # Brown scale: high OM = dark brown, low OM = light tan
        intensity = min(1.0, om_val / 6.0)
        r = 0.76 - 0.45 * intensity
        g = 0.60 - 0.40 * intensity
        b = 0.42 - 0.30 * intensity
        return (r, g, b)

    fig, axes = plt.subplots(1, 3, figsize=(12, max(6, max_depth * 0.06)),
                              gridspec_kw={"width_ratios": [2, 1.5, 1.5]})

    # Panel 1: Profile diagram with texture bars
    ax1 = axes[0]
    for h in profile_data:
        depth = h["bot"] - h["top"]
        color = om_color(h["om"])
        ax1.barh(
            -(h["top"] + depth / 2), 1, height=depth,
            color=color, edgecolor="black", linewidth=0.8
        )
        label = "%s" % h["texture"]
        if h["hzname"]:
            label = "%s: %s" % (h["hzname"], label)
        ax1.text(
            0.5, -(h["top"] + depth / 2), label,
            ha="center", va="center", fontsize=9,
            color="white" if h["om"] > 3 else "black"
        )

    ax1.set_xlim(0, 1)
    ax1.set_ylim(-max_depth, 0)
    ax1.set_ylabel("Depth (cm)")
    ax1.set_title("Soil Profile: %s" % comp_name, fontsize=11, fontweight="bold")
    ax1.set_xticks([])
    ax1.invert_yaxis()
    ax1.set_ylim(-max_depth, 0)

    # Panel 2: Sand/Silt/Clay stacked bars
    ax2 = axes[1]
    for h in profile_data:
        depth = h["bot"] - h["top"]
        y_center = -(h["top"] + depth / 2)
        ax2.barh(y_center, h["sand"], height=depth * 0.8,
                 color="#EDC9AF", edgecolor="black", linewidth=0.3, label="Sand" if h is profile_data[0] else "")
        ax2.barh(y_center, h["silt"], left=h["sand"], height=depth * 0.8,
                 color="#C4A882", edgecolor="black", linewidth=0.3, label="Silt" if h is profile_data[0] else "")
        ax2.barh(y_center, h["clay"], left=h["sand"] + h["silt"], height=depth * 0.8,
                 color="#8B6914", edgecolor="black", linewidth=0.3, label="Clay" if h is profile_data[0] else "")

    ax2.set_xlim(0, 100)
    ax2.set_ylim(-max_depth, 0)
    ax2.set_xlabel("Percent (%)")
    ax2.set_title("Texture (Sand/Silt/Clay)", fontsize=10)
    ax2.legend(loc="lower right", fontsize=8)
    ax2.set_yticks([])

    # Panel 3: pH and OM profiles
    ax3 = axes[2]
    mid_depths = [-(h["top"] + (h["bot"] - h["top"]) / 2) for h in profile_data]
    phs = [h["ph"] for h in profile_data]
    oms = [h["om"] for h in profile_data]

    color_ph = "#1565c0"
    color_om = "#2e7d32"

    ax3.plot(phs, mid_depths, "o-", color=color_ph, linewidth=2, markersize=6, label="pH")
    ax3.set_xlabel("pH", color=color_ph)
    ax3.set_xlim(3, 10)
    ax3.tick_params(axis="x", labelcolor=color_ph)
    ax3.set_title("pH and OM", fontsize=10)
    ax3.set_ylim(-max_depth, 0)
    ax3.set_yticks([])

    ax3b = ax3.twiny()
    ax3b.plot(oms, mid_depths, "s-", color=color_om, linewidth=2, markersize=6, label="OM %")
    ax3b.set_xlabel("Organic Matter (%)", color=color_om)
    ax3b.tick_params(axis="x", labelcolor=color_om)

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "soil_profile.png"), dpi=150, bbox_inches="tight")
    plt.close()

    return True


def generate_report(soil_data, mu_info, lat, lon, mukey, output_dir):
    """Generate markdown report with soil profile table."""
    lines = [
        "# Soil Report",
        "",
        "**Location**: %.4fN, %.4fW  " % (lat, abs(lon)) if lon < 0 else
        "**Location**: %.4fN, %.4fE  " % (lat, lon),
        "**Map Unit Key (mukey)**: %s  " % mukey,
    ]

    if mu_info and len(mu_info) > 0:
        mu = mu_info[0]
        lines.append("**Map Unit Name**: %s  " % mu.get("muname", "N/A"))
        lines.append("**Survey Area**: %s (%s)  " % (
            mu.get("areaname", "N/A"), mu.get("areasymbol", "N/A")))
        lines.append("**Map Unit Kind**: %s  " % mu.get("mukind", "N/A"))
    lines.append("**Generated**: %s  " % datetime.now().strftime("%Y-%m-%d %H:%M"))
    lines.append("")

    if not soil_data:
        lines.append("No soil horizon data returned for this location.")
        report_path = os.path.join(output_dir, "report.md")
        with open(report_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines))
        return

    # Group by component
    components = {}
    for row in soil_data:
        comp = row.get("compname", "Unknown")
        if comp not in components:
            components[comp] = {
                "pct": row.get("comppct_r", ""),
                "taxclass": row.get("taxclname", ""),
                "taxorder": row.get("taxorder", ""),
                "taxsuborder": row.get("taxsuborder", ""),
                "drainage": row.get("drainagecl", ""),
                "horizons": [],
            }
        components[comp]["horizons"].append(row)

    # Component summary
    lines += [
        "## Map Unit Components",
        "",
        "| Component | Pct (%) | Tax. Order | Drainage Class |",
        "|-----------|---------|------------|----------------|",
    ]
    for comp_name, comp_info in components.items():
        lines.append("| %s | %s | %s | %s |" % (
            comp_name,
            safe_float(comp_info["pct"], "%.0f"),
            comp_info["taxorder"] or "--",
            comp_info["drainage"] or "--",
        ))

    # Dominant component details
    dom_name = list(components.keys())[0]
    dom = components[dom_name]
    lines += [
        "",
        "## Dominant Component: %s" % dom_name,
        "",
        "**Taxonomic Classification**: %s  " % (dom["taxclass"] or "N/A"),
        "**Taxonomic Order**: %s  " % (dom["taxorder"] or "N/A"),
        "**Taxonomic Suborder**: %s  " % (dom["taxsuborder"] or "N/A"),
        "**Drainage Class**: %s  " % (dom["drainage"] or "N/A"),
        "",
    ]

    # Drainage interpretation
    drain_interp = interpret_drainage(dom["drainage"])
    lines.append("**Drainage Interpretation**: %s" % drain_interp)
    lines.append("")

    # Horizon table
    lines += [
        "## Soil Profile by Horizon",
        "",
        "| Horizon | Depth (cm) | Sand % | Silt % | Clay % | Texture | OM % | pH | AWC (cm/cm) | Ksat (um/s) |",
        "|---------|-----------|--------|--------|--------|---------|------|-----|-------------|-------------|",
    ]

    for h in dom["horizons"]:
        top = safe_float(h.get("hzdept_r"), "%.0f")
        bot = safe_float(h.get("hzdepb_r"), "%.0f")
        sand = safe_float(h.get("sandtotal_r"))
        silt = safe_float(h.get("silttotal_r"))
        clay = safe_float(h.get("claytotal_r"))
        texture = classify_texture(
            h.get("sandtotal_r"), h.get("silttotal_r"), h.get("claytotal_r")
        )
        om = safe_float(h.get("om_r"))
        ph = safe_float(h.get("ph1to1h2o_r"))
        awc = safe_float(h.get("awc_r"), "%.2f")
        ksat = safe_float(h.get("ksat_r"), "%.2f")
        hz_name = h.get("hzname") or "--"
        depth_str = "%s-%s" % (top, bot)

        lines.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            hz_name, depth_str, sand, silt, clay, texture, om, ph, awc, ksat))

    # Surface horizon interpretations
    surf = dom["horizons"][0] if dom["horizons"] else None
    if surf:
        lines += [
            "",
            "## Surface Horizon Interpretation",
            "",
        ]
        ph_val = surf.get("ph1to1h2o_r")
        om_val = surf.get("om_r")
        lines.append("**pH**: %s -- %s  " % (safe_float(ph_val), interpret_ph(ph_val)))
        lines.append("**Organic Matter**: %s%% -- %s  " % (safe_float(om_val), interpret_om(om_val)))

        awc_val = surf.get("awc_r")
        if awc_val is not None:
            awc_f = float(awc_val)
            if awc_f < 0.10:
                awc_interp = "Low -- frequent irrigation may be needed"
            elif awc_f < 0.18:
                awc_interp = "Moderate -- typical for loamy soils"
            else:
                awc_interp = "High -- good water retention"
            lines.append("**Available Water Capacity**: %s cm/cm -- %s  " % (
                safe_float(awc_val, "%.2f"), awc_interp))

        surf_texture = classify_texture(
            surf.get("sandtotal_r"), surf.get("silttotal_r"), surf.get("claytotal_r")
        )
        lines.append("**Texture**: %s  " % surf_texture)

    # Figures
    lines += [
        "",
        "## Figures",
        "",
        "![Soil Profile](figures/soil_profile.png)",
        "",
    ]

    # References
    lines += [
        "## Data Source and References",
        "",
        "- USDA-NRCS Soil Survey Geographic (SSURGO) Database via Soil Data Access (SDA) API",
        "- Soil Survey Staff. Web Soil Survey. USDA-NRCS. https://websoilsurvey.nrcs.usda.gov/",
        "- Soil Survey Division Staff (2017). Soil Survey Manual. USDA Handbook 18.",
        "- Soil texture classification follows USDA-NRCS texture triangle (Soil Survey Manual, Ch. 3)",
        "",
        "**Note**: SSURGO data represents mapped soil survey estimates, not site-specific samples. "
        "For precision agriculture, supplement with actual soil testing.",
        "",
        "---",
        "*Generated by [AgriClaw claw-soil](https://github.com/bgtamang/AgriClaw) v0.1.0*",
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
        cmd = "python soil.py --demo"
    else:
        cmd = "python soil.py --lat %s --lon %s" % (args_dict["lat"], args_dict["lon"])

    with open(os.path.join(repro_dir, "commands.sh"), "w", encoding="utf-8") as f:
        f.write("#!/bin/bash\n")
        f.write("# Reproduce this analysis\n")
        f.write("# Generated: %s\n\n" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        f.write("pip install matplotlib\n")
        f.write("%s\n" % cmd)

    # environment.yml
    with open(os.path.join(repro_dir, "environment.yml"), "w", encoding="utf-8") as f:
        f.write("name: claw-soil\n")
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
            checksums.append("%s  %s" % (sha256_file(fpath), rel))

    with open(os.path.join(repro_dir, "checksums.sha256"), "w", encoding="utf-8") as f:
        f.write("\n".join(checksums) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="claw-soil: Soil Property Lookup via USDA SDA (AgriClaw)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  python soil.py --demo\n"
               "  python soil.py --lat 40.08 --lon -88.20\n"
               "  python soil.py --lat 42.03 --lon -93.47 --output iowa_soil\n",
    )
    parser.add_argument("--demo", action="store_true",
                        help="Run demo: Champaign IL (40.08N, -88.20W)")
    parser.add_argument("--lat", type=float, help="Latitude (decimal degrees)")
    parser.add_argument("--lon", type=float, help="Longitude (decimal degrees)")
    parser.add_argument("--output", default="output",
                        help="Output directory (default: output)")
    args = parser.parse_args()

    # Demo mode
    if args.demo:
        args.lat = 40.08
        args.lon = -88.20
        args.output = "demo_output"
        print("Running demo: Champaign IL (40.08N, 88.20W)")

    # Validate
    if not args.demo and (args.lat is None or args.lon is None):
        parser.error("Provide --lat/--lon or --demo")

    lat = args.lat
    lon = args.lon

    # Step 1: Get mukey for the point
    print("Querying USDA Soil Data Access for (%.4f, %.4f)..." % (lat, lon))
    mukey = get_mukey_for_point(lat, lon)
    if mukey is None:
        print("ERROR: Could not retrieve map unit key for this location.", file=sys.stderr)
        print("This may mean the location is outside SSURGO coverage (US only).", file=sys.stderr)
        sys.exit(1)
    print("Map unit key (mukey): %s" % mukey)

    # Step 2: Get soil properties
    print("Fetching soil properties...")
    soil_data = get_soil_properties(mukey)
    if soil_data is None:
        print("ERROR: Failed to query soil properties.", file=sys.stderr)
        sys.exit(1)

    if len(soil_data) == 0:
        print("WARNING: No horizon data returned for mukey %s." % mukey, file=sys.stderr)

    # Step 3: Get map unit info
    mu_info = get_mapunit_info(mukey)

    # Create output directory structure
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    tables_dir = os.path.join(output_dir, "tables")
    os.makedirs(tables_dir, exist_ok=True)

    # Write raw horizon data CSV
    if soil_data:
        csv_fields = [
            "compname", "comppct_r", "hzname", "hzdept_r", "hzdepb_r",
            "sandtotal_r", "silttotal_r", "claytotal_r", "om_r",
            "ph1to1h2o_r", "awc_r", "ksat_r", "dbthirdbar_r", "cec7_r",
            "taxclname", "taxorder", "taxsuborder", "drainagecl",
        ]
        # Add texture class column
        enriched = []
        for row in soil_data:
            r = dict(row)
            r["texture_class"] = classify_texture(
                r.get("sandtotal_r"), r.get("silttotal_r"), r.get("claytotal_r")
            )
            enriched.append(r)
        csv_fields.append("texture_class")
        write_csv(os.path.join(tables_dir, "soil_horizons.csv"), enriched, csv_fields)
        print("Wrote soil_horizons.csv (%d rows)" % len(enriched))

    # Write component summary CSV
    if soil_data:
        comp_seen = {}
        comp_rows = []
        for row in soil_data:
            cn = row.get("compname", "")
            if cn not in comp_seen:
                comp_seen[cn] = True
                comp_rows.append({
                    "compname": cn,
                    "comppct_r": row.get("comppct_r", ""),
                    "taxclname": row.get("taxclname", ""),
                    "taxorder": row.get("taxorder", ""),
                    "drainagecl": row.get("drainagecl", ""),
                })
        write_csv(
            os.path.join(tables_dir, "component_summary.csv"),
            comp_rows,
            ["compname", "comppct_r", "taxclname", "taxorder", "drainagecl"],
        )
        print("Wrote component_summary.csv (%d components)" % len(comp_rows))

    # Generate profile plot
    plots_ok = generate_profile_plot(soil_data, output_dir)
    if plots_ok:
        print("Generated figure: soil_profile.png")

    # Generate report
    generate_report(soil_data, mu_info, lat, lon, mukey, output_dir)
    print("Generated report.md")

    # Reproducibility bundle
    args_dict = {"demo": args.demo, "lat": lat, "lon": lon}
    generate_repro_bundle(output_dir, args_dict)
    print("Generated reproducibility bundle (commands.sh, environment.yml, checksums.sha256)")

    # Summary
    print("")
    print("=" * 50)
    if mu_info and len(mu_info) > 0:
        print("Map Unit: %s" % mu_info[0].get("muname", "N/A"))
        print("Survey Area: %s" % mu_info[0].get("areaname", "N/A"))

    if soil_data:
        dom = soil_data[0]
        print("Dominant Component: %s (%s%%)" % (
            dom.get("compname", "N/A"), safe_float(dom.get("comppct_r"), "%.0f")))
        print("Drainage: %s" % (dom.get("drainagecl") or "N/A"))
        print("Taxonomic Order: %s" % (dom.get("taxorder") or "N/A"))
        surf_texture = classify_texture(
            dom.get("sandtotal_r"), dom.get("silttotal_r"), dom.get("claytotal_r")
        )
        print("Surface Texture: %s (Sand: %s%%, Silt: %s%%, Clay: %s%%)" % (
            surf_texture,
            safe_float(dom.get("sandtotal_r")),
            safe_float(dom.get("silttotal_r")),
            safe_float(dom.get("claytotal_r")),
        ))
        print("Surface pH: %s" % safe_float(dom.get("ph1to1h2o_r")))
        print("Surface OM: %s%%" % safe_float(dom.get("om_r")))
        print("Horizons: %d" % sum(
            1 for h in soil_data if h.get("compname") == dom.get("compname")))
    else:
        print("No soil data returned for this location.")

    print("")
    print("Output: %s/" % os.path.abspath(output_dir))


if __name__ == "__main__":
    main()
