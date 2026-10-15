#!/usr/bin/env python3
"""claw-soil: mapped soil properties for a US location from USDA-NRCS SSURGO via the
Soil Data Access (SDA) REST API: texture (SSURGO's own class), component-weighted
map-unit properties, available water storage by depth, drainage and water table,
plus a USDA texture-triangle classifier for sand/silt/clay data.

Usage:
    python scripts/soil.py --demo
    python scripts/soil.py --lat 40.08 --lon -88.20
    python scripts/soil.py --mukey 242963
    python scripts/soil.py --texture-input lab_samples.csv

Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT licence.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import platform
import re
import shlex
import sys
import urllib.error
import urllib.request
from datetime import datetime
from pathlib import Path

SKILL = "claw-soil"
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
# USDA texture classes (Soil Science Division Staff 2017, Soil Survey Manual, pp. 122-124)
# ---------------------------------------------------------------------------

TEXTURE_CLASSES = ["sand", "loamy sand", "sandy loam", "loam", "silt loam", "silt", "sandy clay loam",
                   "clay loam", "silty clay loam", "sandy clay", "silty clay", "clay"]
SUM_TOLERANCE = 2.0  # sand + silt + clay may differ from 100 by this much (rounding in source data)


def _texture_rules(sand: float, silt: float, clay: float) -> dict[str, bool]:
    """The 12 class definitions, word for word from the Soil Survey Manual (2017)."""
    s15 = round(silt + 1.5 * clay, 9)
    s2 = round(silt + 2.0 * clay, 9)
    return {
        "sand": sand > 85 and s15 < 15,
        "loamy sand": 70 <= sand <= 90 and s15 >= 15 and s2 < 30,
        "sandy loam": (7 <= clay < 20 and sand > 52 and s2 >= 30) or (clay < 7 and silt < 50 and s2 >= 30),
        "loam": 7 <= clay < 27 and 28 <= silt < 50 and sand <= 52,
        "silt loam": (silt >= 50 and 12 <= clay < 27) or (50 <= silt < 80 and clay < 12),
        "silt": silt >= 80 and clay < 12,
        "sandy clay loam": 20 <= clay < 35 and silt < 28 and sand > 45,
        "clay loam": 27 <= clay < 40 and 20 < sand <= 45,
        "silty clay loam": 27 <= clay < 40 and sand <= 20,
        "sandy clay": clay >= 35 and sand > 45,
        "silty clay": clay >= 40 and silt >= 40,
        "clay": clay >= 40 and sand <= 45 and silt < 40,
    }


def usda_texture_class(sand, silt, clay, tol: float = SUM_TOLERANCE) -> tuple[str | None, str]:
    """USDA texture class from % sand, silt and clay of the fine-earth fraction.

    Returns (class or None, note). Sums within 100 +/- tol are rescaled to 100 first.
    """
    try:
        vals = [float(v) for v in (sand, silt, clay)]
    except (TypeError, ValueError):
        return None, "sand, silt or clay missing"
    if any(v < 0 or v > 100 for v in vals):
        return None, "a fraction is outside 0-100%"
    total = sum(vals)
    if abs(total - 100.0) > tol:
        return None, f"sand + silt + clay = {total:g}%, not 100 ± {tol:g}"
    sand, silt, clay = (round(v * 100.0 / total, 9) for v in vals)
    hits = [k for k, ok in _texture_rules(sand, silt, clay).items() if ok]
    if len(hits) != 1:  # cannot happen if the definitions partition the triangle (tested)
        return None, f"ambiguous: {hits}"
    note = "" if abs(total - 100.0) < 1e-9 else f"rescaled from a sum of {total:g}%"
    return hits[0], note


def base_texture_from_texdesc(texdesc: str | None) -> str | None:
    """Reduce an SSURGO texture description to one of the 12 basic classes.

    'Gravelly silt loam' -> 'silt loam'; 'Loamy fine sand' -> 'loamy sand'.
    Stratified or 'x to y' descriptions and non-mineral materials return None.
    """
    if not texdesc:
        return None
    s = texdesc.strip().lower()
    if " to " in s or "stratified" in s:
        return None
    s = re.sub(r"\b(?:very\s+)?(?:coarse|fine)\s+(?=sand)", "", s)
    s = re.sub(r"\s+", " ", s).strip()
    for cls in sorted(TEXTURE_CLASSES, key=len, reverse=True):
        if s == cls or s.endswith(" " + cls):
            return cls
    return None


# ---------------------------------------------------------------------------
# Interpretation tables
# ---------------------------------------------------------------------------

def reaction_class(ph) -> str:
    """Soil Survey Manual (2017) Table 3-15 reaction classes (pH rounded to 0.1)."""
    if ph is None:
        return ""
    p = round(float(ph), 1)
    for upper, name in ((3.4, "Ultra acid"), (4.4, "Extremely acid"), (5.0, "Very strongly acid"),
                        (5.5, "Strongly acid"), (6.0, "Moderately acid"), (6.5, "Slightly acid"),
                        (7.3, "Neutral"), (7.8, "Slightly alkaline"), (8.4, "Moderately alkaline"),
                        (9.0, "Strongly alkaline")):
        if p <= upper:
            return name
    return "Very strongly alkaline"


# Soil Survey Manual (2017), pp. 213-215, shortened.
DRAINAGE_NOTES = {
    "excessively drained": "Water is removed very rapidly; commonly coarse textured or very shallow.",
    "somewhat excessively drained": "Water is removed rapidly; commonly coarse textured or very shallow.",
    "well drained": "Water is removed readily but not rapidly; wetness does not inhibit root growth for "
                    "significant periods in most growing seasons.",
    "moderately well drained": "Wet for only a short time within the rooting depth during the growing "
                               "season, but long enough that most mesophytic crops are affected.",
    "somewhat poorly drained": "Wet at a shallow depth for significant periods; wetness markedly restricts "
                               "mesophytic crops unless artificially drained.",
    "poorly drained": "Wet at shallow depths periodically or for long periods; most mesophytic crops cannot "
                      "be grown unless the soil is artificially drained.",
    "very poorly drained": "Free water at or very near the surface during much of the growing season; most "
                           "mesophytic crops cannot be grown unless artificially drained.",
    "subaqueous": "Free water is above the soil surface.",
}

# Organic matter guidance: adjustable rule of thumb (see SKILL.md), not a sourced classification.
OM_GUIDE = [(1.0, "very low"), (2.0, "low"), (4.0, "moderate"), (6.0, "high")]


def om_guidance(om) -> str:
    if om is None:
        return ""
    for upper, name in OM_GUIDE:
        if float(om) < upper:
            return name
    return "very high"


# ---------------------------------------------------------------------------
# Soil Data Access
# ---------------------------------------------------------------------------

SDA_URL = "https://sdmdataaccess.sc.egov.usda.gov/tabular/post.rest"
SQL_POINT = "SELECT * FROM SDA_Get_Mukey_from_intersection_with_WktWgs84('point({lon:.6f} {lat:.6f})')"
SQL_MAPUNIT = """SELECT mu.mukey, mu.musym, mu.muname, mu.mukind, l.areasymbol, l.areaname, sc.saverest,
 ma.aws025wta, ma.aws050wta, ma.aws0100wta, ma.aws0150wta, ma.drclassdcd, ma.hydgrpdcd,
 ma.wtdepannmin, ma.wtdepaprjunmin, ma.brockdepmin, ma.slopegradwta, ma.niccdcd, ma.flodfreqdcd
FROM mapunit mu
INNER JOIN legend l ON l.lkey = mu.lkey
LEFT JOIN sacatalog sc ON sc.areasymbol = l.areasymbol
LEFT JOIN muaggatt ma ON ma.mukey = mu.mukey
WHERE mu.mukey = '{mukey}'"""
SQL_HORIZONS = """SELECT co.cokey, co.compname, co.comppct_r, co.majcompflag, co.compkind, co.taxclname, co.taxorder,
 co.drainagecl, co.hydgrp, ch.chkey, ch.hzname, ch.hzdept_r, ch.hzdepb_r, ch.sandtotal_r, ch.silttotal_r,
 ch.claytotal_r, ch.om_r, ch.ph1to1h2o_r, ch.awc_r, ch.ksat_r, ch.dbthirdbar_r, ch.cec7_r, ch.wthirdbar_r,
 ch.wfifteenbar_r, tg.texdesc
FROM component co
LEFT JOIN chorizon ch ON ch.cokey = co.cokey
LEFT JOIN chtexturegrp tg ON tg.chkey = ch.chkey AND tg.rvindicator = 'Yes'
WHERE co.mukey = '{mukey}'
ORDER BY co.comppct_r DESC, co.cokey, ch.hzdept_r"""

DEMO = {"lat": 40.08, "lon": -88.20, "point": "sda_point_40.08_-88.20.json",
        "mapunit": "sda_mapunit_242963.json", "horizons": "sda_horizons_242963.json"}


def sda_query(sql: str, timeout: int = 120) -> dict:
    payload = json.dumps({"query": sql, "format": "JSON+COLUMNNAME"}).encode("utf-8")
    req = urllib.request.Request(SDA_URL, data=payload, headers={
        "Content-Type": "application/json", "User-Agent": f"AgriClaw/{VERSION} ({SKILL})"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            text = resp.read().decode("utf-8")
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"Soil Data Access returned HTTP {e.code}: "
                           f"{e.read().decode('utf-8', 'replace')[:300]}") from None
    except (urllib.error.URLError, TimeoutError) as e:
        raise RuntimeError(f"Could not reach USDA Soil Data Access ({e}).") from None
    try:
        return json.loads(text) if text.strip() else {}
    except json.JSONDecodeError:
        raise RuntimeError(f"Soil Data Access returned non-JSON text: {text[:300]}") from None


def parse_sda_table(resp: dict) -> list[dict]:
    """SDA 'JSON+COLUMNNAME' -> list of dicts. An empty result is {} or a header-only table."""
    if not resp:
        return []
    table = resp.get("Table")
    if table is None:
        raise RuntimeError(f"Unexpected Soil Data Access response: {str(resp)[:300]}")
    if len(table) < 2:
        return []
    return [dict(zip(table[0], row)) for row in table[1:]]


def check_mukey(mukey) -> str:
    m = str(mukey).strip()
    if not m.isdigit():
        raise ValueError(f"mukey '{mukey}' must be digits only.")
    return m


def num(v):
    if v is None or (isinstance(v, str) and not v.strip()):
        return None
    return float(v)


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------


def build_components(rows: list[dict]) -> list[dict]:
    """Group SDA horizon rows into components (ordered by comppct_r, as returned)."""
    comps: dict[str, dict] = {}
    for r in rows:
        c = comps.setdefault(str(r["cokey"]), {
            "cokey": str(r["cokey"]), "compname": r["compname"], "comppct_r": num(r["comppct_r"]) or 0.0,
            "majcompflag": (r.get("majcompflag") or "").strip(), "compkind": r.get("compkind") or "",
            "taxclname": r.get("taxclname") or "", "taxorder": r.get("taxorder") or "",
            "drainagecl": r.get("drainagecl") or "", "hydgrp": r.get("hydgrp") or "", "horizons": []})
        if r.get("chkey") is None or r.get("hzdept_r") is None or r.get("hzdepb_r") is None:
            continue
        h = {"hzname": r.get("hzname") or "", "top_cm": num(r["hzdept_r"]), "bottom_cm": num(r["hzdepb_r"]),
             "sand_pct": num(r.get("sandtotal_r")), "silt_pct": num(r.get("silttotal_r")),
             "clay_pct": num(r.get("claytotal_r")), "om_pct": num(r.get("om_r")),
             "ph_h2o": num(r.get("ph1to1h2o_r")), "awc_cm_cm": num(r.get("awc_r")),
             "ksat_um_s": num(r.get("ksat_r")), "db_g_cm3": num(r.get("dbthirdbar_r")),
             "cec7_meq_100g": num(r.get("cec7_r")), "w33kpa_pct_vol": num(r.get("wthirdbar_r")),
             "w1500kpa_pct_vol": num(r.get("wfifteenbar_r")), "texdesc": r.get("texdesc") or ""}
        cls, note = usda_texture_class(h["sand_pct"], h["silt_pct"], h["clay_pct"])
        h["texture_computed"] = cls or ""
        h["texture_computed_note"] = note
        base = base_texture_from_texdesc(h["texdesc"])
        h["texture_class"] = h["texdesc"] or (cls.capitalize() if cls else "")
        h["texture_source"] = ("SSURGO texdesc" if h["texdesc"] else
                               "computed from sand/silt/clay" if cls else "not available")
        h["texture_agrees"] = ("" if not (base and cls) else "yes" if base == cls else "no")
        h["ph_class"] = reaction_class(h["ph_h2o"])
        c["horizons"].append(h)
    return list(comps.values())


def depth_weighted(horizons: list[dict], col: str, top: float, bottom: float) -> tuple[float | None, float]:
    """Thickness-weighted mean of col over [top, bottom) cm, and the thickness with data."""
    tot = thick = 0.0
    for h in horizons:
        ov = min(h["bottom_cm"], bottom) - max(h["top_cm"], top)
        if ov > 0 and h[col] is not None:
            tot += h[col] * ov
            thick += ov
    return (tot / thick if thick else None), thick


def mapunit_weighted(components: list[dict], col: str, top: float, bottom: float):
    """Component-percent-weighted mean of the depth-weighted values (components with data only)."""
    num_, den = 0.0, 0.0
    for c in components:
        v, thick = depth_weighted(c["horizons"], col, top, bottom)
        if v is not None:
            num_ += c["comppct_r"] * v
            den += c["comppct_r"]
    return (num_ / den if den else None), den


def aws_cm(horizons: list[dict], depth_cm: float) -> float:
    """Available water storage (cm) from the surface to depth_cm: sum of AWC x thickness."""
    return sum(h["awc_cm_cm"] * max(0.0, min(h["bottom_cm"], depth_cm) - max(h["top_cm"], 0.0))
               for h in horizons if h["awc_cm_cm"] is not None)


def aws_profile(components: list[dict], depths=range(10, 210, 10)) -> list[dict]:
    """Component-weighted cumulative AWS (mm) by depth, the same aggregation SSURGO uses for
    muaggatt.aws*wta (weighted average of components)."""
    with_hz = [c for c in components if c["horizons"]]
    den = sum(c["comppct_r"] for c in with_hz)
    dom = with_hz[0] if with_hz else None
    out = []
    for d in depths:
        mu = sum(c["comppct_r"] * aws_cm(c["horizons"], d) for c in with_hz) / den if den else None
        out.append({"depth_cm": d, "aws_mm": round(mu * 10.0, 2) if mu is not None else "",
                    "aws_mm_dominant_component": round(aws_cm(dom["horizons"], d) * 10.0, 2) if dom else "",
                    "max_described_depth_cm": max((h["bottom_cm"] for c in with_hz for h in c["horizons"]),
                                                  default="")})
    return out


WEIGHTED_PROPS = [("sand_pct", "%"), ("silt_pct", "%"), ("clay_pct", "%"), ("om_pct", "%"),
                  ("ph_h2o", "pH (1:1 water)"), ("awc_cm_cm", "cm/cm"), ("db_g_cm3", "g/cm3"),
                  ("cec7_meq_100g", "meq/100 g"), ("w33kpa_pct_vol", "% vol"), ("w1500kpa_pct_vol", "% vol")]


def weighted_table(components: list[dict], ranges=((0, 30), (0, 100))) -> list[dict]:
    out = []
    for top, bottom in ranges:
        for col, unit in WEIGHTED_PROPS:
            v, pct = mapunit_weighted(components, col, top, bottom)
            out.append({"property": col, "depth_range_cm": f"{top}-{bottom}", "value": "" if v is None else round(v, 3),
                        "units": unit, "pct_of_mapunit_with_data": round(pct, 1)})
    return out


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def make_figure(components: list[dict], profile: list[dict], mu: dict, out: Path) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dom = next((c for c in components if c["horizons"]), None)
    if dom is None:
        return []
    hz = dom["horizons"]
    depth = max(h["bottom_cm"] for h in hz)
    fig, axes = plt.subplots(1, 3, figsize=(12, 6), gridspec_kw={"width_ratios": [2.2, 1.2, 1.4]})
    ax = axes[0]
    for h in hz:
        th = h["bottom_cm"] - h["top_cm"]
        left = 0.0
        for col, color in (("sand_pct", "#e3c79b"), ("silt_pct", "#b79a6b"), ("clay_pct", "#7a5a2e")):
            v = h[col] or 0.0
            ax.barh(h["top_cm"] + th / 2, v, height=th * 0.95, left=left, color=color, edgecolor="white", lw=0.5)
            left += v
        ax.text(101, h["top_cm"] + th / 2, f"{h['hzname']}: {h['texture_class'] or 'n/a'}", va="center", fontsize=7)
    ax.set_xlim(0, 175)
    ax.set_ylim(depth, 0)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("sand | silt | clay (%)")
    ax.set_ylabel("Depth (cm)")
    ax.set_title(f"{dom['compname']} ({dom['comppct_r']:.0f}% of map unit)", fontsize=10)
    ax = axes[1]
    mids = [(h["top_cm"] + h["bottom_cm"]) / 2 for h in hz]
    ax.plot([h["ph_h2o"] if h["ph_h2o"] is not None else float("nan") for h in hz], mids, "o-", color="#2e6da4")
    ax.set_xlabel("pH (1:1 water)", color="#2e6da4")
    ax.set_ylim(depth, 0)
    ax2 = ax.twiny()
    ax2.plot([h["om_pct"] if h["om_pct"] is not None else float("nan") for h in hz], mids, "s--", color="#2f6b3a")
    ax2.set_xlabel("Organic matter (%)", color="#2f6b3a")
    ax = axes[2]
    pd_ = [p for p in profile if p["aws_mm"] != ""]
    ax.plot([p["aws_mm"] for p in pd_], [p["depth_cm"] for p in pd_], "-", color="#1f5f8b", lw=2,
            label="computed (all components)")
    pts = [(25, mu.get("aws025wta")), (50, mu.get("aws050wta")), (100, mu.get("aws0100wta")),
           (150, mu.get("aws0150wta"))]
    pts = [(d, num(v) * 10) for d, v in pts if num(v) is not None]
    if pts:
        ax.plot([p[1] for p in pts], [p[0] for p in pts], "o", color="#c0392b", label="SSURGO muaggatt")
    ax.set_ylim(max(depth, 150), 0)
    ax.set_xlabel("Cumulative available water (mm)")
    ax.legend(fontsize=7, loc="lower left")
    ax.grid(alpha=0.3)
    fig.suptitle(mu.get("muname", ""), fontsize=11)
    fig.tight_layout()
    (out / "figures").mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "figures" / "soil_profile.png", dpi=150)
    plt.close(fig)
    return ["figures/soil_profile.png"]


REFERENCES = [
    "Soil Science Division Staff (2017). Soil survey manual. C. Ditzler, K. Scheffe & H.C. Monger (eds.). "
    "USDA Handbook 18. Government Printing Office, Washington, D.C. Texture classes pp. 122-124; reaction "
    "classes Table 3-15; drainage classes pp. 213-215. "
    "https://www.nrcs.usda.gov/resources/guides-and-instructions/soil-survey-manual",
    "Soil Survey Staff, Natural Resources Conservation Service, USDA. Soil Survey Geographic (SSURGO) Database, "
    "accessed through Soil Data Access. https://sdmdataaccess.sc.egov.usda.gov/",
    "USDA-NRCS (2025). Soil Data Mart 4.3 Table Column Descriptions (muaggatt, chorizon, chtexturegrp). "
    "https://sdmdataaccess.nrcs.usda.gov/documents/TableColumnDescriptionsReport.pdf",
    "Allen, R.G., Pereira, L.S., Raes, D. & Smith, M. (1998). Crop evapotranspiration: guidelines for "
    "computing crop water requirements. FAO Irrigation and Drainage Paper 56. FAO, Rome (TAW, Chapter 8). "
    "https://www.fao.org/4/x0490e/x0490e00.htm",
]


def write_report(out: Path, ctx: dict) -> None:
    mu = ctx["mapunit"]
    L = [f"# {SKILL} report", ""]
    if ctx["warnings"]:
        L += ["## Warnings", ""] + [f"- {w}" for w in ctx["warnings"]] + [""]
    L += ["## Map unit", "",
          f"- Location: {ctx['location']}",
          f"- Map unit: {mu.get('muname', '')} (symbol {mu.get('musym', '')}, mukey {mu.get('mukey', '')}, "
          f"{mu.get('mukind', '')})",
          f"- Survey area: {mu.get('areaname', '')} ({mu.get('areasymbol', '')}), version date {mu.get('saverest', '')}",
          "", "SSURGO is mapped soil-survey data: each map unit is a mix of components drawn at roughly "
          "1:12,000-1:24,000 scale, and property values are representative estimates for each component. It is "
          "not a sample of this field. Use soil tests and on-site observation for decisions.", "",
          "## Components", "", "| Component | % of map unit | Major | Drainage class | Hydrologic group | Taxonomic class |",
          "|---|---|---|---|---|---|"]
    for c in ctx["components"]:
        L.append(f"| {c['compname']} | {c['comppct_r']:.0f} | {c['majcompflag']} | {c['drainagecl'] or 'n/a'} | "
                 f"{c['hydgrp'] or 'n/a'} | {c['taxclname'] or 'n/a'} |")
    dom = next((c for c in ctx["components"] if c["horizons"]), None)
    if dom:
        L += ["", f"## Horizons of the dominant component ({dom['compname']})", "",
              "Texture is SSURGO's own representative texture class (`chtexturegrp.texdesc`, RV row). The "
              "'computed' column applies the USDA texture triangle to sand/silt/clay as a cross-check.", "",
              "| Horizon | Depth (cm) | Texture (SSURGO) | Computed class | Sand/Silt/Clay (%) | OM (%) | pH (class) | AWC (cm/cm) |",
              "|---|---|---|---|---|---|---|---|"]
        for h in dom["horizons"]:
            ssc = "/".join("n/a" if h[k] is None else f"{h[k]:g}" for k in ("sand_pct", "silt_pct", "clay_pct"))
            L.append(f"| {h['hzname']} | {h['top_cm']:g}-{h['bottom_cm']:g} | {h['texdesc'] or 'n/a'} | "
                     f"{h['texture_computed'] or 'n/a'} | {ssc} | {'n/a' if h['om_pct'] is None else h['om_pct']} | "
                     f"{'n/a' if h['ph_h2o'] is None else h['ph_h2o']} ({h['ph_class'] or 'n/a'}) | "
                     f"{'n/a' if h['awc_cm_cm'] is None else h['awc_cm_cm']} |")
    w = ctx["water"]
    L += ["", "## Water storage and drainage (for claw-water-balance)", "",
          "| Item | Value | Source |", "|---|---|---|",
          f"| Available water storage 0-25 / 0-50 / 0-100 / 0-150 cm | {w['aws_0_25cm_mm']} / {w['aws_0_50cm_mm']} / "
          f"{w['aws_0_100cm_mm']} / {w['aws_0_150cm_mm']} mm | SSURGO muaggatt (weighted average of components) |",
          f"| Same, recomputed from horizon AWC | {w['aws_0_25cm_mm_computed']} / {w['aws_0_50cm_mm_computed']} / "
          f"{w['aws_0_100cm_mm_computed']} / {w['aws_0_150cm_mm_computed']} mm | this skill (check) |",
          f"| Drainage class (dominant condition) | {w['drainage_class']} | muaggatt.drclassdcd |",
          f"| Shallowest water table, any month / April-June | {w['water_table_min_cm']} / {w['water_table_apr_jun_min_cm']} cm | "
          "muaggatt (components ≥ 15%) |",
          f"| Shallowest bedrock | {str(w['bedrock_min_cm']) + ' cm' if w['bedrock_min_cm'] else 'none recorded'} | "
          "muaggatt.brockdepmin |",
          f"| Hydrologic soil group | {w['hydrologic_group']} | muaggatt.hydgrpdcd |", "",
          "Total available water for a root depth Zr: interpolate `tables/aws_profile.csv` at Zr (FAO-56 TAW = "
          "1000 (θFC − θWP) Zr; SSURGO AWC is θFC − θWP adjusted for rock fragments and salinity).", ""]
    if dom and dom["drainagecl"]:
        note = DRAINAGE_NOTES.get(dom["drainagecl"].lower(), "")
        L += ["## General guidance (not a recommendation)", "",
              f"- Drainage, {dom['drainagecl']}: {note} (Soil Survey Manual 2017 definition.)"]
        surf = dom["horizons"][0] if dom["horizons"] else None
        if surf and surf["ph_h2o"] is not None:
            L.append(f"- Surface pH {surf['ph_h2o']}: {surf['ph_class']} (Soil Survey Manual Table 3-15). Lime or "
                     "acidification needs depend on the crop and a current soil test, not on mapped values.")
        if surf and surf["om_pct"] is not None:
            L.append(f"- Surface organic matter {surf['om_pct']}%: {om_guidance(surf['om_pct'])} by an adjustable "
                     "rule of thumb (<1 very low, 1-2 low, 2-4 moderate, 4-6 high, ≥6 very high).")
        L.append("")
    L += ["## Figures", ""] + [f"![{f}]({f})" for f in ctx["figures"]] + [""]
    L += ["## References", ""] + [f"- {r}" for r in REFERENCES] + ["", f"*{SKILL} {VERSION}*"]
    write_text(out / "report.md", "\n".join(L) + "\n")


def run_texture_csv(path: str, out: Path) -> tuple[list[dict], list[str]]:
    p = Path(path)
    if not p.is_file():
        raise ValueError(f"Input file not found: {path}")
    with p.open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError("The texture CSV has no rows.")
    need = [c for c in ("sand_pct", "silt_pct", "clay_pct") if c not in rows[0]]
    if need:
        raise ValueError(f"Missing column(s) {', '.join(need)}; found {', '.join(rows[0].keys())}. Use % of the "
                         "fine-earth (<2 mm) fraction.")
    out_rows, warnings = [], []
    for r in rows:
        cls, note = usda_texture_class(r["sand_pct"] or None, r["silt_pct"] or None, r["clay_pct"] or None)
        out_rows.append(dict(r, usda_texture_class=cls or "", note=note))
    bad = sum(1 for r in out_rows if not r["usda_texture_class"])
    if bad:
        warnings.append(f"{bad} row(s) could not be classified (see the 'note' column).")
    return out_rows, warnings


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="soil.py", description=__doc__.splitlines()[0])
    ap.add_argument("--demo", action="store_true", help="Drummer silty clay loam near Champaign, IL (bundled SDA data)")
    ap.add_argument("--lat", type=float, help="latitude, decimal degrees (WGS84)")
    ap.add_argument("--lon", type=float, help="longitude, decimal degrees (WGS84)")
    ap.add_argument("--mukey", help="SSURGO map unit key (instead of --lat/--lon)")
    ap.add_argument("--texture-input", help="CSV with sand_pct, silt_pct, clay_pct to classify (no SSURGO lookup)")
    ap.add_argument("--out", default=f"agriclaw_output/{SKILL}", help="output folder")
    return ap


def main(argv: list[str] | None = None, query=sda_query) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = build_parser()
    args = ap.parse_args(argv)
    out = Path(args.out)
    warnings: list[str] = []
    try:
        if args.texture_input:
            out.mkdir(parents=True, exist_ok=True)
            rows, warnings = run_texture_csv(args.texture_input, out)
            save_table(rows, out, "texture_classes.csv")
            counts: dict[str, int] = {}
            for r in rows:
                counts[r["usda_texture_class"] or "unclassified"] = counts.get(r["usda_texture_class"] or "unclassified", 0) + 1
            L = [f"# {SKILL} report: texture classes", ""]
            if warnings:
                L += ["## Warnings", ""] + [f"- {w}" for w in warnings] + [""]
            L += ["USDA texture classes from sand/silt/clay (Soil Survey Manual 2017, pp. 122-124). Sums within "
                  f"100 ± {SUM_TOLERANCE:g}% were rescaled to 100.", "", "| Class | Samples |", "|---|---|"]
            L += [f"| {k} | {v} |" for k, v in sorted(counts.items())]
            L += ["", "## References", "", f"- {REFERENCES[0]}", "", f"*{SKILL} {VERSION}*"]
            write_text(out / "report.md", "\n".join(L) + "\n")
            write_repro_bundle(out, argv, inputs=[args.texture_input], packages=["matplotlib"])
            print(f"Done: {len(rows)} samples classified. Report: {(out / 'report.md').resolve()}")
            return 0

        raw = {}
        inputs = []
        if args.demo:
            def load(name):
                p = SKILL_DIR / "demo" / name
                inputs.append(str(p))
                return json.loads(p.read_text(encoding="utf-8"))
            pt = parse_sda_table(load(DEMO["point"]))
            mukey = check_mukey(pt[0]["mukey"])
            mu_rows = parse_sda_table(load(DEMO["mapunit"]))
            hz_rows = parse_sda_table(load(DEMO["horizons"]))
            location = f"{DEMO['lat']}, {DEMO['lon']} (bundled SDA responses)"
        else:
            if args.mukey:
                mukey = check_mukey(args.mukey)
                location = f"mukey {mukey}"
            elif args.lat is not None and args.lon is not None:
                if not (-90 <= args.lat <= 90 and -180 <= args.lon <= 180):
                    ap.error("--lat must be within -90..90 and --lon within -180..180")
                raw["point"] = query(SQL_POINT.format(lat=args.lat, lon=args.lon))
                pt = parse_sda_table(raw["point"])
                if not pt:
                    raise ValueError(f"No SSURGO map unit at {args.lat}, {args.lon}. SSURGO covers the USA and some "
                                     "territories; check that latitude/longitude are not swapped and the point is on land.")
                mukey = check_mukey(pt[0]["mukey"])
                if len(pt) > 1:
                    warnings.append(f"The point touches {len(pt)} map units; the first ({mukey}) is used.")
                location = f"{args.lat}, {args.lon} (WGS84)"
            else:
                ap.error("give --lat and --lon, --mukey, --texture-input or --demo")
            raw["mapunit"] = query(SQL_MAPUNIT.format(mukey=mukey))
            raw["horizons"] = query(SQL_HORIZONS.format(mukey=mukey))
            mu_rows, hz_rows = parse_sda_table(raw["mapunit"]), parse_sda_table(raw["horizons"])
            write_text(out / "raw" / "sda_responses.json", json.dumps(raw, indent=1))
        if not mu_rows:
            raise ValueError(f"mukey {mukey} was not found in SSURGO.")
        mu = mu_rows[0]
        components = build_components(hz_rows)
        if not components or not any(c["horizons"] for c in components):
            warnings.append("No horizon data for this map unit (e.g. water, urban land or rock outcrop).")
        pct_total = sum(c["comppct_r"] for c in components)
        no_hz = [c for c in components if not c["horizons"]]
        if no_hz:
            warnings.append(f"Components without horizon data ({', '.join(c['compname'] for c in no_hz)}; "
                            f"{sum(c['comppct_r'] for c in no_hz):.0f}% of the map unit) are left out of the "
                            "weighted values.")
        if components and abs(pct_total - 100) > 1:
            warnings.append(f"Component percentages sum to {pct_total:g}%, not 100%.")
        disagree = [(c["compname"], h["hzname"], h["texdesc"], h["texture_computed"]) for c in components
                    for h in c["horizons"] if h["texture_agrees"] == "no"]
        if disagree:
            warnings.append(f"{len(disagree)} horizon(s) where SSURGO's texture class differs from the class computed "
                            "from its sand/silt/clay (e.g. {} {}: '{}' vs '{}'); SSURGO's class is reported."
                            .format(*disagree[0]))
        profile = aws_profile(components)
        water = {"mukey": mukey}
        for d, key in ((25, "aws025wta"), (50, "aws050wta"), (100, "aws0100wta"), (150, "aws0150wta")):
            v = num(mu.get(key))
            water[f"aws_0_{d}cm_mm"] = round(v * 10.0, 1) if v is not None else ""
            comp = next((p["aws_mm"] for p in aws_profile(components, [d])), "")
            water[f"aws_0_{d}cm_mm_computed"] = round(comp, 1) if comp != "" else ""
            if v is not None and comp != "" and abs(comp - v * 10) > max(2.0, 0.05 * v * 10):
                warnings.append(f"Recomputed AWS 0-{d} cm ({comp:.1f} mm) differs from SSURGO muaggatt "
                                f"({v * 10:.1f} mm); check the horizon data.")
        water.update({
            "drainage_class": mu.get("drclassdcd") or "", "hydrologic_group": mu.get("hydgrpdcd") or "",
            "water_table_min_cm": mu.get("wtdepannmin") or "", "water_table_apr_jun_min_cm": mu.get("wtdepaprjunmin") or "",
            "bedrock_min_cm": mu.get("brockdepmin") or "", "slope_pct": mu.get("slopegradwta") or "",
            "nonirrigated_capability_class": mu.get("niccdcd") or "", "flooding": mu.get("flodfreqdcd") or "",
            "muname": mu.get("muname") or "", "areasymbol": mu.get("areasymbol") or "",
            "survey_version_date": mu.get("saverest") or ""})
        wt = num(mu.get("wtdepannmin"))
        if wt is not None and wt < 200:
            warnings.append(f"A water table can be as shallow as {wt:g} cm (muaggatt.wtdepannmin). FAO-56 assumes "
                            "capillary rise is about zero only when the water table is more than ~1 m below the "
                            "root zone; tell claw-water-balance users, and ask whether the field is tile-drained.")
        bedrock = num(mu.get("brockdepmin"))
        if bedrock is not None and bedrock < 150:
            warnings.append(f"Bedrock as shallow as {bedrock:g} cm limits rooting depth and water storage.")

        hz_out = [{"compname": c["compname"], "comppct_r": c["comppct_r"], "majcompflag": c["majcompflag"],
                   **{k: v for k, v in h.items()}} for c in components for h in c["horizons"]]
        save_table(hz_out or [{"compname": ""}], out, "horizons.csv")
        save_table([{k: v for k, v in c.items() if k != "horizons"} | {"n_horizons": len(c["horizons"])}
                    for c in components] or [{"compname": ""}], out, "components.csv")
        save_table(weighted_table(components), out, "mapunit_weighted_properties.csv")
        save_table(profile, out, "aws_profile.csv")
        save_table([water], out, "mapunit_summary.csv")
        figures = make_figure(components, profile, mu, out)
        write_report(out, {"warnings": warnings, "mapunit": mu, "components": components, "water": water,
                           "location": location, "figures": figures})
        write_repro_bundle(out, argv, inputs=inputs, packages=["matplotlib"])
    except (ValueError, RuntimeError) as e:
        print(ascii_safe(f"Error: {e}"), file=sys.stderr)
        return 1
    print(ascii_safe(f"Done: {mu.get('muname')} (mukey {mukey}); AWS 0-100 cm {water['aws_0_100cm_mm']} mm. "
                     f"Report: {(out / 'report.md').resolve()}"))
    for w in warnings:
        print(ascii_safe(f"Warning: {w}"), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
