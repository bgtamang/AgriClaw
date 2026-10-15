"""Tests for claw-soil. All offline: Soil Data Access responses are saved in tests/fixtures/."""

import csv
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "soil.py"
FIX = SKILL_DIR / "tests" / "fixtures"
spec = importlib.util.spec_from_file_location("claw_soil", SCRIPT)
soil = importlib.util.module_from_spec(spec)
spec.loader.exec_module(soil)


def fixture(name):
    return json.loads((FIX / name).read_text(encoding="utf-8"))


def cls(sand, silt, clay):
    return soil.usda_texture_class(sand, silt, clay)[0]


# --- USDA texture triangle (Soil Survey Manual 2017) ----------------------------------------

@pytest.mark.parametrize("sand,silt,clay,expected", [
    (92, 5, 3, "sand"), (82, 12, 6, "loamy sand"), (65, 25, 10, "sandy loam"), (40, 40, 20, "loam"),
    (20, 65, 15, "silt loam"), (5, 88, 7, "silt"), (60, 13, 27, "sandy clay loam"), (33, 34, 33, "clay loam"),
    (10, 56, 34, "silty clay loam"), (52, 6, 42, "sandy clay"), (6, 47, 47, "silty clay"), (20, 20, 60, "clay"),
])
def test_one_point_in_each_of_the_12_classes(sand, silt, clay, expected):
    assert cls(sand, silt, clay) == expected


def test_points_the_old_classifier_got_wrong():
    assert cls(10, 70, 20) == "silt loam"      # old code: silty clay loam
    assert cls(60, 35, 5) == "sandy loam"      # old code: loamy sand


@pytest.mark.parametrize("sand,silt,clay,expected", [
    (52, 30, 18, "loam"),               # sand <= 52 is loam; sandy loam needs > 52
    (53, 29, 18, "sandy loam"),
    (30, 43, 27, "clay loam"),          # clay 27 leaves loam (< 27)
    (20, 40, 40, "silty clay"),         # clay 40 and silt 40 -> silty clay, not clay (silt < 40)
    (21, 39, 40, "clay"),
    (9, 80, 11, "silt"),                # silt >= 80 with clay < 12
    (10, 79, 11, "silt loam"),
    (45, 20, 35, "clay loam"),          # sand 45 is not > 45
    (46, 19, 35, "sandy clay"),
    (20, 50, 30, "silty clay loam"),    # sand <= 20
    (21, 49, 30, "clay loam"),
    (90, 0, 10, "loamy sand"),          # silt + 1.5 clay = 15 is >= 15
    (86, 14, 0, "sand"),                # silt + 1.5 clay = 14 < 15
    (85, 15, 0, "loamy sand"),          # sand 85 is not > 85
    (70, 30, 0, "sandy loam"),          # silt + 2 clay = 30 is >= 30
    (71, 29, 0, "loamy sand"),
    (52, 41, 7, "loam"),                # clay 7 is in 7..<27
    (28, 60, 12, "silt loam"),          # clay 12 with silt >= 50
    (30, 50, 20, "silt loam"),          # silt 50 is >= 50
    (52, 28, 20, "loam"),               # silt 28 is not < 28 (sandy clay loam)
    (53, 27, 20, "sandy clay loam"),
])
def test_class_boundaries(sand, silt, clay, expected):
    assert cls(sand, silt, clay) == expected


def test_definitions_partition_the_triangle():
    """Every grid point (0.5% steps) satisfies exactly one of the 12 definitions."""
    n = 0
    for i in range(0, 201):
        for j in range(0, 201 - i):
            sand, clay = i / 2, j / 2
            silt = 100 - sand - clay
            hits = [k for k, ok in soil._texture_rules(sand, silt, clay).items() if ok]
            assert len(hits) == 1, (sand, silt, clay, hits)
            n += 1
    assert n == 20301


def test_sum_tolerance():
    c, note = soil.usda_texture_class(30, 40, 29)            # sum 99: rescaled, then classified
    assert c == "clay loam" and "rescaled" in note
    c, note = soil.usda_texture_class(30, 40, 20)            # sum 90: refused
    assert c is None and "not 100" in note
    assert soil.usda_texture_class(None, 40, 20)[0] is None


@pytest.mark.parametrize("texdesc,expected", [
    ("Silty clay loam", "silty clay loam"), ("Gravelly silt loam", "silt loam"), ("Loamy fine sand", "loamy sand"),
    ("Fine sandy loam", "sandy loam"), ("Very fine sandy loam", "sandy loam"), ("Coarse sand", "sand"),
    ("Mucky loamy sand", "loamy sand"), ("Loam", "loam"), ("Clay", "clay"), ("Silty clay", "silty clay"),
    ("Stratified sandy loam to clay loam", None), ("Highly decomposed plant material", None),
    ("Unweathered bedrock", None), ("", None),
])
def test_base_class_from_ssurgo_texdesc(texdesc, expected):
    assert soil.base_texture_from_texdesc(texdesc) == expected


@pytest.mark.parametrize("ph,expected", [
    (3.4, "Ultra acid"), (3.5, "Extremely acid"), (4.4, "Extremely acid"), (4.5, "Very strongly acid"),
    (5.0, "Very strongly acid"), (5.1, "Strongly acid"), (6.0, "Moderately acid"), (6.5, "Slightly acid"),
    (6.6, "Neutral"), (7.3, "Neutral"), (7.4, "Slightly alkaline"), (8.4, "Moderately alkaline"),
    (9.0, "Strongly alkaline"), (9.1, "Very strongly alkaline"),
])
def test_reaction_classes_table_3_15(ph, expected):
    assert soil.reaction_class(ph) == expected


# --- SDA parsing and aggregation ------------------------------------------------------------

def test_parse_point_and_empty_responses():
    assert soil.parse_sda_table(fixture("sda_point_40.08_-88.20.json")) == [{"mukey": "242963"}]
    assert soil.parse_sda_table(fixture("sda_point_outside_coverage.json")) == []
    with pytest.raises(RuntimeError):
        soil.parse_sda_table({"unexpected": 1})
    with pytest.raises(ValueError):
        soil.check_mukey("242963'; DROP TABLE mapunit; --")


@pytest.mark.parametrize("mukey", ["242963", "426520"])
def test_aws_reproduces_ssurgo_muaggatt(mukey):
    """Component-weighted AWS recomputed from horizon AWC equals SSURGO's muaggatt values."""
    mu = soil.parse_sda_table(fixture(f"sda_mapunit_{mukey}.json"))[0]
    comps = soil.build_components(soil.parse_sda_table(fixture(f"sda_horizons_{mukey}.json")))
    for d, key in ((25, "aws025wta"), (50, "aws050wta"), (100, "aws0100wta"), (150, "aws0150wta")):
        computed = soil.aws_profile(comps, [d])[0]["aws_mm"]
        assert computed == pytest.approx(float(mu[key]) * 10.0, abs=0.1), (mukey, d)


def test_ssurgo_texture_is_primary_and_agrees():
    comps = soil.build_components(soil.parse_sda_table(fixture("sda_horizons_426520.json")))
    hz = [h for c in comps for h in c["horizons"]]
    organic = [h for h in hz if h["hzname"] in ("Oe", "Oa")]
    assert organic and all(h["texture_source"] == "SSURGO texdesc" and h["texture_computed"] == "" for h in organic)
    mucky = next(h for h in hz if h["texdesc"] == "Mucky loamy sand")
    assert mucky["texture_computed"] == "loamy sand" and mucky["texture_agrees"] == "yes"
    assert all(h["texture_agrees"] in ("yes", "") for h in hz)
    drummer = soil.build_components(soil.parse_sda_table(fixture("sda_horizons_242963.json")))
    assert [c["compname"] for c in drummer] == ["Drummer", "Peotone", "Harpster"]
    assert drummer[0]["horizons"][0]["texture_class"] == "Silty clay loam"


def test_component_weighting_hand_calc():
    def comp(pct, clay, top=0, bottom=100):
        return {"comppct_r": pct, "horizons": [{"top_cm": top, "bottom_cm": bottom, "clay_pct": clay}]}
    comps = [comp(60, 20.0), comp(30, 40.0), comp(10, None)]
    v, pct = soil.mapunit_weighted(comps, "clay_pct", 0, 30)
    assert v == pytest.approx((60 * 20 + 30 * 40) / 90) and pct == 90
    two_layer = [{"top_cm": 0, "bottom_cm": 20, "clay_pct": 10.0}, {"top_cm": 20, "bottom_cm": 60, "clay_pct": 30.0}]
    assert soil.depth_weighted(two_layer, "clay_pct", 0, 30)[0] == pytest.approx((10 * 20 + 30 * 10) / 30)


# --- Command line ---------------------------------------------------------------------------

def fake_query(mukey, point_file):
    def q(sql):
        if "SDA_Get_Mukey" in sql:
            return fixture(point_file)
        if "FROM mapunit" in sql:
            return fixture(f"sda_mapunit_{mukey}.json")
        return fixture(f"sda_horizons_{mukey}.json")
    return q


def test_point_outside_coverage_fails_with_advice(tmp_path):
    rc = soil.main(["--lat", "42.0", "--lon", "-87.0", "--out", str(tmp_path)],
                   query=fake_query("242963", "sda_point_outside_coverage.json"))
    assert rc == 1


def test_lookup_by_mukey_with_saved_responses(tmp_path):
    rc = soil.main(["--mukey", "426520", "--out", str(tmp_path)], query=fake_query("426520", "unused"))
    assert rc == 0
    with open(tmp_path / "tables" / "mapunit_summary.csv", newline="", encoding="utf-8") as f:
        s = next(csv.DictReader(f))
    assert s["muname"].startswith("Leola loamy sand") and float(s["aws_0_100cm_mm"]) == pytest.approx(98.8)
    assert (tmp_path / "raw" / "sda_responses.json").is_file()


def test_texture_csv_mode(tmp_path):
    with open(tmp_path / "lab.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["sample", "sand_pct", "silt_pct", "clay_pct"])
        w.writerows([["a", 10, 70, 20], ["b", 60, 35, 5], ["c", 30, 30, 20]])
    assert soil.main(["--texture-input", str(tmp_path / "lab.csv"), "--out", str(tmp_path / "o")]) == 0
    with open(tmp_path / "o" / "tables" / "texture_classes.csv", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert [r["usda_texture_class"] for r in rows] == ["silt loam", "sandy loam", ""]
    assert "not 100" in rows[2]["note"]


def test_cli_demo(tmp_path):
    r = subprocess.run([sys.executable, str(SCRIPT), "--demo", "--out", str(tmp_path)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    for f in ["report.md", "tables/horizons.csv", "tables/components.csv", "tables/mapunit_weighted_properties.csv",
              "tables/aws_profile.csv", "tables/mapunit_summary.csv", "figures/soil_profile.png",
              "reproducibility/commands.sh", "reproducibility/environment.txt", "reproducibility/checksums.sha256"]:
        assert (tmp_path / f).is_file(), f
    rep = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "Drummer silty clay loam" in rep and "not a sample of this field" in rep
    assert rep.index("## Warnings") < rep.index("## Map unit") and "water table" in rep
    with open(tmp_path / "tables" / "aws_profile.csv", newline="", encoding="utf-8") as f:
        prof = {int(r["depth_cm"]): float(r["aws_mm"]) for r in csv.DictReader(f)}
    assert prof[100] == pytest.approx(173.9, abs=0.1)
