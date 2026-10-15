"""Tests for claw-yield-units: hand-checkable conversions, published worked examples,
round trips, plot-to-area scaling, batch mode and the CLI demo."""

import csv
import importlib.util
import random
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "yield_units.py"
spec = importlib.util.spec_from_file_location("claw_yield_units", SCRIPT)
yu = importlib.util.module_from_spec(spec)
spec.loader.exec_module(yu)


# --- unit definitions and bushel weights ------------------------------------------------------

def test_exact_unit_definitions():
    assert yu.ACRE_M2 == pytest.approx(4046.8564224, abs=1e-9)
    assert yu.LB_AC_TO_KG_HA == pytest.approx(1.12085116, abs=1e-8)
    assert yu.LENGTH_TO_M["in"] * 12 == pytest.approx(yu.LENGTH_TO_M["ft"])


def test_bushel_weights_match_ah697_table5():
    """AH-697 Table 5: 60-lb bushel = 27.2155 kg; 56-lb = 0.0254 t; 48-lb = 0.021772 t; 32-lb = 0.014515 t."""
    assert yu.CROPS["soybean"]["lb_per_bu"] * yu.LB_KG == pytest.approx(27.2155, abs=1e-4)
    assert yu.CROPS["wheat"]["lb_per_bu"] == 60 and yu.CROPS["corn"]["lb_per_bu"] == 56
    assert yu.CROPS["corn"]["lb_per_bu"] * yu.LB_KG / 1000 == pytest.approx(0.0254, abs=1e-5)
    assert yu.CROPS["barley"]["lb_per_bu"] * yu.LB_KG / 1000 == pytest.approx(0.021772, abs=1e-6)
    assert yu.CROPS["oats"]["lb_per_bu"] * yu.LB_KG / 1000 == pytest.approx(0.014515, abs=1e-6)
    for c in ("sorghum", "rye", "flaxseed"):
        assert yu.CROPS[c]["lb_per_bu"] == 56
    assert yu.CROPS["rice"]["lb_per_bu"] == 45


def test_standard_moistures():
    assert yu.CROPS["corn"]["std_moisture_pct"] == 15.5
    assert yu.CROPS["soybean"]["std_moisture_pct"] == 13.0
    assert yu.CROPS["wheat"]["std_moisture_pct"] == 13.5
    assert yu.CROPS["sorghum"]["std_moisture_pct"] is None        # sources disagree: ask the user


# --- hand-checkable conversions -----------------------------------------------------------------

def test_soybean_60_bu_ac():
    # 60 bu/ac x 60 lb/bu = 3600 lb/ac = 4035.06 kg/ha
    assert yu.convert(60, "bu/ac", "lb/ac", 60) == pytest.approx(3600.0)
    assert yu.convert(60, "bu/ac", "kg/ha", 60) == pytest.approx(4035.064, abs=0.001)
    assert yu.convert(60, "bu/ac", "t/ha", 60) == pytest.approx(4.035064, abs=1e-6)


def test_corn_200_bu_ac():
    # 200 bu/ac x 56 lb/bu = 11200 lb/ac = 12553.53 kg/ha
    assert yu.convert(200, "bu/ac", "kg/ha", 56) == pytest.approx(12553.53, abs=0.01)
    assert yu.convert(12.0, "t/ha", "bu/ac", 56) == pytest.approx(191.18, abs=0.01)


def test_other_units():
    assert yu.convert(1, "t/ha", "g/m2") == pytest.approx(100.0)
    assert yu.convert(1, "cwt/ac", "lb/ac") == pytest.approx(100.0)
    assert yu.convert(1000, "kg/ha", "lb/ac") == pytest.approx(892.179, abs=0.001)
    with pytest.raises(ValueError):
        yu.convert(1, "bu/ac", "kg/ha")                 # no bushel weight
    with pytest.raises(ValueError):
        yu.convert(1, "bushels", "kg/ha")


def test_round_trips():
    rng = random.Random(1)
    units = ["kg/ha", "t/ha", "g/m2", "lb/ac", "cwt/ac", "bu/ac"]
    for _ in range(300):
        v = rng.uniform(0.1, 20000)
        a, b = rng.choice(units), rng.choice(units)
        back = yu.convert(yu.convert(v, a, b, 56), b, a, 56)
        assert back == pytest.approx(v, rel=1e-12)


# --- moisture adjustment ------------------------------------------------------------------------

def test_uf_ag442_example():
    """UF/IFAS AG442: 10,000 lb/ac corn at 20% moisture = 9,467 lb/ac at 15.5%."""
    assert yu.adjust_moisture(10000, 20, 15.5) == pytest.approx(9467.46, abs=0.01)


def test_arkansas_fsa1078_examples():
    """FSA-1078: 261,000 lb of soybean dry matter = 300,000 lb at 13%, 310,714 lb at 16%, 290,000 lb at 10%."""
    assert yu.adjust_moisture(261000, 0, 13) == pytest.approx(300000, abs=0.5)
    assert yu.adjust_moisture(261000, 0, 16) == pytest.approx(310714.3, abs=0.5)
    assert yu.adjust_moisture(261000, 0, 10) == pytest.approx(290000, abs=0.5)
    assert yu.adjust_moisture(300000, 13, 16) == pytest.approx(310714.3, abs=0.5)


def test_moisture_identities_and_errors():
    assert yu.adjust_moisture(5000, 15, 15) == 5000
    assert yu.adjust_moisture(yu.adjust_moisture(5000, 22, 15.5), 15.5, 22) == pytest.approx(5000)
    with pytest.raises(ValueError):
        yu.adjust_moisture(5000, 100, 15)
    with pytest.raises(ValueError):
        yu.adjust_moisture(5000, -1, 15)


# --- plot to area -------------------------------------------------------------------------------

def test_plot_area_imperial():
    # 2 rows x 30 in x 20 ft = 100 ft2 = 9.290304 m2; 10 lb -> 4356 lb/ac -> 72.6 bu/ac soybean
    a = yu.plot_area_m2(2, 30, 20, "in", "ft")
    assert a == pytest.approx(9.290304, abs=1e-9)
    kg_ha = yu.plot_yield_kg_ha(10, "lb", a)
    assert kg_ha / yu.LB_AC_TO_KG_HA == pytest.approx(4356.0, abs=1e-6)
    assert yu.from_kg_ha(kg_ha, "bu/ac", 60) == pytest.approx(72.6, abs=1e-9)


def test_plot_area_metric_and_direct_area():
    a = yu.plot_area_m2(2, 0.76, 5.0, "m", "m")
    assert a == pytest.approx(7.6)
    assert yu.plot_yield_kg_ha(2.5, "kg", a) == pytest.approx(3289.47, abs=0.01)
    assert yu.plot_area_m2(area=100, area_unit="ft2") == pytest.approx(9.290304)
    assert yu.plot_yield_kg_ha(2500, "g", 7.6) == pytest.approx(3289.47, abs=0.01)
    with pytest.raises(ValueError):
        yu.plot_area_m2(2, 30, None, "in", "ft")


# --- command line -------------------------------------------------------------------------------

def read_yields(folder):
    with open(Path(folder) / "tables" / "yields.csv", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def test_cli_single_value_with_moisture(tmp_path):
    assert yu.main(["--value", "10000", "--from", "lb/ac", "--crop", "corn", "--moisture", "20",
                    "--out", str(tmp_path)]) == 0
    r = read_yields(tmp_path)[0]
    assert float(r["yield_lb_ac"]) == pytest.approx(9467.456, abs=0.001)
    assert r["moisture_basis"] == "15.5% moisture" and r["lb_per_bu"] == "56.0"
    rep = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "56 lb/bu (USDA AH-697)" in rep and "15.5% (UF/IFAS AG442" in rep


def test_cli_plot_requires_units(tmp_path):
    assert yu.main(["--plot-weight", "10", "--weight-unit", "lb", "--rows", "2", "--row-spacing", "30",
                    "--length", "20", "--crop", "soybean", "--out", str(tmp_path)]) == 1
    assert yu.main(["--plot-weight", "10", "--weight-unit", "lb", "--rows", "2", "--row-spacing", "30",
                    "--row-spacing-unit", "in", "--length", "20", "--length-unit", "ft", "--crop", "soybean",
                    "--moisture", "13", "--out", str(tmp_path)]) == 0
    assert float(read_yields(tmp_path)[0]["yield_bu_ac"]) == pytest.approx(72.6, abs=0.001)


def test_cli_sorghum_needs_target_moisture(tmp_path):
    base = ["--value", "8000", "--from", "kg/ha", "--crop", "sorghum", "--moisture", "18", "--out", str(tmp_path)]
    assert yu.main(base) == 1
    assert yu.main(base + ["--to-moisture", "14"]) == 0
    assert float(read_yields(tmp_path)[0]["yield_kg_ha"]) == pytest.approx(8000 * 82 / 86, abs=0.001)


def test_cli_batch_per_plot_geometry_and_moisture(tmp_path):
    with open(tmp_path / "p.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["plot", "plot_weight", "harvest_moisture_pct", "rows", "row_spacing", "length"])
        w.writerow([1, 10, 13, 2, 30, 20])
        w.writerow([2, 10, 18, 4, 15, 20])
    assert yu.main(["--input", str(tmp_path / "p.csv"), "--crop", "soybean", "--weight-unit", "lb",
                    "--row-spacing-unit", "in", "--length-unit", "ft", "--out", str(tmp_path / "o")]) == 0
    rows = read_yields(tmp_path / "o")
    assert float(rows[0]["yield_bu_ac"]) == pytest.approx(72.6, abs=0.001)
    assert float(rows[1]["yield_bu_ac"]) == pytest.approx(72.6 * 82 / 87, abs=0.001)
    assert rows[1]["plot"] == "2" and float(rows[1]["area_m2"]) == pytest.approx(9.2903, abs=1e-4)


def test_cli_flags_implausible_geometry(tmp_path):
    assert yu.main(["--plot-weight", "4", "--weight-unit", "kg", "--rows", "2", "--row-spacing", "30",
                    "--row-spacing-unit", "m", "--length", "5", "--length-unit", "m", "--crop", "corn",
                    "--moisture", "18", "--out", str(tmp_path)]) == 0
    assert "row spacing" in read_yields(tmp_path)[0]["flags"]


def test_cli_demo(tmp_path):
    r = subprocess.run([sys.executable, str(SCRIPT), "--demo", "--out", str(tmp_path)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    for f in ["report.md", "tables/yields.csv", "figures/yields.png", "reproducibility/commands.sh",
              "reproducibility/environment.txt", "reproducibility/checksums.sha256"]:
        assert (tmp_path / f).is_file(), f
    rows = read_yields(tmp_path)
    assert len(rows) == 12
    area_ft2 = 2 * 2.5 * 17.5
    for r in rows:     # independent recomputation of every plot
        lb_ac = float(r["plot_weight"]) * 43560 / area_ft2 * (100 - float(r["harvest_moisture_pct"])) / 87
        assert float(r["yield_bu_ac"]) == pytest.approx(lb_ac / 60, abs=0.001)
    rep = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "60 lb/bu" in rep and "13% (" in rep
