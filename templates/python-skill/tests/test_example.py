"""Tests for claw-example (template for AgriClaw Python tests): known answers,
a bad-input check, and an end-to-end run of the command line on the demo data."""

import csv
import importlib.util
import math
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "example.py"
spec = importlib.util.spec_from_file_location("claw_example", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_known_answers():
    rows = [{"genotype": g, "yield_kg_ha": str(y)}
            for g, y in [("A", 1), ("A", 2), ("A", 3), ("B", 10), ("B", 20)]]
    res = mod.summarise_trait(rows, "yield_kg_ha", "genotype", min_n=3)
    assert [r["mean"] for r in res] == [2, 15]
    assert math.isclose(res[0]["sd"], 1.0)
    assert math.isclose(res[1]["sd"], math.sqrt(50))
    assert [r["low_n"] for r in res] == [False, True]


def test_missing_column_is_named():
    with pytest.raises(ValueError, match="biomass_g"):
        mod.summarise_trait([{"genotype": "A", "yield_kg_ha": "1"}], "biomass_g", "genotype")


def test_cli_demo(tmp_path):
    r = subprocess.run([sys.executable, str(SCRIPT), "--demo", "--out", str(tmp_path)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    for f in ["report.md", "tables/summary.csv", "figures/means.png",
              "reproducibility/commands.sh", "reproducibility/environment.txt",
              "reproducibility/checksums.sha256"]:
        assert (tmp_path / f).is_file(), f
    with (tmp_path / "tables" / "summary.csv").open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 4 and all(r["n"] == "3" for r in rows)
