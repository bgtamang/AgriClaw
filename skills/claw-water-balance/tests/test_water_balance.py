"""Tests for claw-water-balance: FAO-56 worked examples, exact mass balance, stress onset,
irrigation scheduling and an end-to-end demo run. All offline."""

import csv
import importlib.util
import random
import subprocess
import sys
from datetime import date, timedelta
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "water_balance.py"
spec = importlib.util.spec_from_file_location("claw_water_balance", SCRIPT)
wb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wb)


def days(n, kc=1.0, et0=5.0, precip=0.0, start=date(2025, 6, 1), **extra):
    return [dict({"date": start + timedelta(days=i), "kc": kc, "et0_mm": et0, "precip_mm": precip}, **extra)
            for i in range(n)]


# --- FAO-56 worked examples -----------------------------------------------------------------

def test_fao56_example37_water_stress():
    """FAO-56 Example 37: tomato, Zr 0.8 m, p 0.40, thetaFC 0.32, thetaWP 0.12, Dr0 55 mm,
    ET0 5 mm/d, Kc 1.2, no rain: TAW 160, RAW 64 and the published daily Ks and Dr."""
    taw = wb.taw_from_theta(0.32, 0.12, 0.8)
    assert taw == pytest.approx(160.0) and 0.40 * taw == pytest.approx(64.0)
    rows = wb.run_balance(days(10, kc=1.2, et0=5.0), taw, 0.40, dr0=55.0, p_adjust=False)
    ks_pub = [1.00, 1.00, 0.97, 0.91, 0.85, 0.80, 0.75, 0.70, 0.66, 0.62]
    dr_pub = [61.0, 67.0, 72.8, 78.3, 83.4, 88.2, 92.6, 96.9, 100.8, 104.5]
    etc_adj_pub = [6.0, 6.0, 5.8, 5.4, 5.1, 4.8, 4.5, 4.2, 3.9, 3.7]
    for r, ks, dr, e in zip(rows, ks_pub, dr_pub, etc_adj_pub):
        assert r["ks"] == pytest.approx(ks, abs=0.006)
        assert r["dr_end_mm"] == pytest.approx(dr, abs=0.06)
        assert r["etc_adj_mm"] == pytest.approx(e, abs=0.06)


@pytest.mark.parametrize("theta_fc,theta_wp,zr,p,taw,raw", [
    (0.15, 0.06, 1.2, 0.55, 108, 59), (0.32, 0.15, 1.2, 0.55, 204, 112), (0.35, 0.23, 1.2, 0.55, 144, 79),
    (0.15, 0.06, 0.4, 0.30, 36, 11), (0.32, 0.15, 0.8, 0.40, 136, 54),
])
def test_fao56_example36_taw_raw(theta_fc, theta_wp, zr, p, taw, raw):
    """FAO-56 Example 36 (maize, onion and tomato on loamy sand, silt and silty clay)."""
    t = wb.taw_from_theta(theta_fc, theta_wp, zr)
    assert t == pytest.approx(taw, abs=0.5) and p * t == pytest.approx(raw, abs=0.5)


def test_fao56_example27_and_box15_kc_climate_adjustment():
    # Example 27: maize Kc_mid 1.20, h 2 m -> 1.07 at Taipei (u2 1.3, RHmin 75), 1.30 at Mocha (4.6, 44)
    assert wb.kc_climate_adjust(1.20, 1.3, 75, 2.0) == pytest.approx(1.07, abs=0.005)
    assert wb.kc_climate_adjust(1.20, 4.6, 44, 2.0) == pytest.approx(1.30, abs=0.005)
    # Box 15: dry bean Kc_mid 1.15, h 0.4 m, RHmin 30%, u2 2.2 m/s -> 1.19
    assert wb.kc_climate_adjust(1.15, 2.2, 30, 0.4) == pytest.approx(1.19, abs=0.005)
    # standard climate leaves the table value unchanged; limits on u2 and RHmin are applied
    assert wb.kc_climate_adjust(1.15, 2.0, 45, 1.0) == pytest.approx(1.15)
    assert wb.kc_climate_adjust(1.15, 9.0, 5, 1.0) == pytest.approx(wb.kc_climate_adjust(1.15, 6.0, 20, 1.0))


def test_fao56_example28_kc_curve():
    """FAO-56 Example 28: dry bean, stages 25/25/30/20, Kc 0.15/1.19/0.35."""
    st = (25, 25, 30, 20)
    expected = {20: 0.15, 40: 0.77, 70: 1.19, 95: 0.56}
    for i, kc in expected.items():
        assert wb.kc_on_day(i, st, 0.15, 1.19, 0.35)[0] == pytest.approx(kc, abs=0.005)
    assert wb.kc_on_day(40, st, 0.15, 1.19, 0.35)[1] == "development"
    assert wb.kc_on_day(101, st, 0.15, 1.19, 0.35) == (None, "")
    assert wb.kc_on_day(100, st, 0.15, 1.19, 0.35)[0] == pytest.approx(0.35)


def test_p_adjustment_for_etc():
    assert wb.p_adjusted(0.55, 5.0) == pytest.approx(0.55)
    assert wb.p_adjusted(0.55, 8.0) == pytest.approx(0.43)
    assert wb.p_adjusted(0.55, 2.0) == pytest.approx(0.67)
    assert wb.p_adjusted(0.55, 20.0) == 0.1          # lower limit
    assert wb.p_adjusted(0.70, 0.0) == 0.8           # upper limit


# --- Water balance behaviour ----------------------------------------------------------------

def test_mass_balance_closes_exactly():
    rng = random.Random(2025)
    ds = []
    for i in range(200):
        p = rng.choice([0, 0, 0, 0, rng.uniform(0, 5), rng.uniform(5, 60)])
        ds.append({"date": date(2025, 4, 1) + timedelta(days=i), "kc": rng.uniform(0.2, 1.25),
                   "et0_mm": rng.uniform(0.5, 9), "precip_mm": p, "runoff_mm": p * rng.uniform(0, 0.3),
                   "capillary_rise_mm": rng.choice([0, 0, 0.5]),
                   "irrigation_mm": rng.choice([0] * 30 + [25])})
    for kwargs in ({}, {"auto_irrigate": True}, {"auto_irrigate": True, "mad": 0.3}, {"p_adjust": False},
                   {"min_rain_frac_et0": 0.2}):
        rows = wb.run_balance(ds, 150.0, 0.55, dr0=40.0, **kwargs)
        assert abs(wb.mass_balance_residual(rows, 40.0)) < 1e-9
        assert all(-1e-9 <= r["dr_end_mm"] <= 150.0 + 1e-9 for r in rows)
        assert all(0 <= r["ks"] <= 1 for r in rows)


def test_stress_onset_and_decline_analytic():
    """No rain, ETc 6 mm/d, TAW 100, p 0.5 (RAW 50): Dr = 6, 12, ..., 54; stress starts on the
    first day that begins with Dr > RAW (day 10), then TAW - Dr shrinks by 1 - 6/50 per day."""
    rows = wb.run_balance(days(20, kc=1.0, et0=6.0), 100.0, 0.5, p_adjust=False)
    assert [r["stressed"] for r in rows[:9]] == [False] * 9 and rows[9]["stressed"]
    assert rows[8]["dr_end_mm"] == pytest.approx(54.0)
    for k in range(9, 20):
        assert 100 - rows[k]["dr_end_mm"] == pytest.approx(46.0 * 0.88 ** (k - 8), rel=1e-9)
    assert rows[9]["ks"] == pytest.approx(0.92)


def test_rain_refills_and_drains():
    ds = days(3, kc=1.0, et0=5.0)
    ds[1]["precip_mm"] = 30.0
    rows = wb.run_balance(ds, 100.0, 0.5, dr0=10.0, p_adjust=False)
    # day 1: Dr 10 -> 15; day 2: 30 mm rain refills 15 mm, 15 mm drains, then ET 5 -> Dr 5
    assert rows[0]["dr_end_mm"] == pytest.approx(15.0)
    assert rows[1]["deep_percolation_mm"] == pytest.approx(15.0) and rows[1]["dr_end_mm"] == pytest.approx(5.0)
    assert rows[2]["dr_end_mm"] == pytest.approx(10.0)


def test_runoff_and_ignored_small_rain():
    ds = days(2, kc=1.0, et0=5.0)
    ds[0].update(precip_mm=20.0, runoff_mm=5.0)
    ds[1].update(precip_mm=0.8)
    rows = wb.run_balance(ds, 100.0, 0.5, dr0=50.0, p_adjust=False, min_rain_frac_et0=0.2)
    assert rows[0]["eff_precip_mm"] == pytest.approx(15.0) and rows[0]["dr_end_mm"] == pytest.approx(40.0)
    assert rows[1]["precip_ignored_mm"] == pytest.approx(0.8) and rows[1]["dr_end_mm"] == pytest.approx(45.0)


def test_auto_irrigation_at_raw_prevents_stress():
    rows = wb.run_balance(days(90, kc=1.1, et0=6.0), 120.0, 0.55, auto_irrigate=True)
    events = [r for r in rows if r["auto_irrigation_mm"] > 0]
    assert events and not any(r["stressed"] for r in rows)
    for r in events:
        assert r["auto_irrigation_mm"] == pytest.approx(r["dr_start_mm"])     # refill to field capacity
        assert r["dr_start_mm"] >= r["raw_mm"] - 1e-9
        assert r["deep_percolation_mm"] == pytest.approx(0.0)


def test_et_limited_when_storage_runs_out():
    rows = wb.run_balance(days(3, kc=1.0, et0=15.0), 20.0, 0.5, p_adjust=False)
    assert rows[1]["et_limited_by_storage"] and rows[1]["dr_end_mm"] == pytest.approx(20.0)
    assert rows[2]["etc_adj_mm"] == pytest.approx(0.0)
    assert abs(wb.mass_balance_residual(rows, 0.0)) < 1e-12


def test_soil_profile_interpolation(tmp_path):
    (tmp_path / "tables").mkdir()
    with open(tmp_path / "tables" / "aws_profile.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["depth_cm", "aws_mm", "aws_mm_dominant_component", "max_described_depth_cm"])
        for d in range(10, 210, 10):
            w.writerow([d, 1.5 * d, 1.5 * d, 152])
    prof, summ = wb.read_soil_dir(str(tmp_path))
    assert wb.interp_aws(prof, 95) == pytest.approx(142.5) and wb.interp_aws(prof, 5) == pytest.approx(7.5)
    assert summ["max_described_depth_cm"] == 152.0
    with pytest.raises(ValueError):
        wb.interp_aws(prof, 250)


# --- Command line ---------------------------------------------------------------------------

def write_weather(path, n=160, start=date(2025, 5, 1), hole=None):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["date", "precip_mm", "et0_mm"])
        for i in range(n):
            d = start + timedelta(days=i)
            w.writerow([d.isoformat(), "" if d == hole else ("12" if i % 9 == 0 else "0"), "4.5"])


def test_cli_requires_one_soil_source(tmp_path):
    write_weather(tmp_path / "w.csv")
    base = ["--weather", str(tmp_path / "w.csv"), "--crop", "soybean", "--plant-date", "2025-05-10",
            "--out", str(tmp_path / "o")]
    assert wb.main(base) == 1                                         # no soil information: refuse
    assert wb.main(base + ["--taw-mm", "150", "--theta-fc", "0.3", "--theta-wp", "0.1"]) == 1
    assert wb.main(base + ["--taw-mm", "150"]) == 0


def test_cli_missing_weather_day_fails(tmp_path):
    write_weather(tmp_path / "w.csv", hole=date(2025, 6, 15))
    rc = wb.main(["--weather", str(tmp_path / "w.csv"), "--crop", "soybean", "--plant-date", "2025-05-10",
                  "--taw-mm", "150", "--out", str(tmp_path / "o")])
    assert rc == 1


def test_cli_irrigated_run_writes_events(tmp_path):
    write_weather(tmp_path / "w.csv")
    rc = wb.main(["--weather", str(tmp_path / "w.csv"), "--crop", "custom", "--kc-ini", "0.4", "--kc-mid", "1.15",
                  "--kc-end", "0.5", "--stage-days", "20,30,60,25", "--root-depth-m", "0.8", "--p", "0.5",
                  "--theta-fc", "0.30", "--theta-wp", "0.14", "--plant-date", "2025-05-10", "--irrigate",
                  "--out", str(tmp_path / "o")])
    assert rc == 0
    with open(tmp_path / "o" / "tables" / "season_summary.csv", newline="", encoding="utf-8") as f:
        s = next(csv.DictReader(f))
    assert int(s["irrigation_events"]) > 0 and s["stress_days"] == "0"
    assert float(s["taw_mm"]) == pytest.approx(1000 * 0.16 * 0.8)


def test_cli_demo(tmp_path):
    r = subprocess.run([sys.executable, str(SCRIPT), "--demo", "--out", str(tmp_path)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    for f in ["report.md", "tables/daily_water_balance.csv", "tables/season_summary.csv", "tables/crop_parameters.csv",
              "tables/irrigation_events.csv", "figures/depletion.png", "figures/kc_et.png",
              "reproducibility/commands.sh", "reproducibility/environment.txt", "reproducibility/checksums.sha256"]:
        assert (tmp_path / f).is_file(), f
    with open(tmp_path / "tables" / "season_summary.csv", newline="", encoding="utf-8") as f:
        s = next(csv.DictReader(f))
    assert s["days"] == "150" and abs(float(s["mass_balance_residual_mm"])) < 1e-6
    # Kc_mid recomputed independently from the demo weather (mid-season = days 71-120 after planting)
    with open(SKILL_DIR / "demo" / "daily_weather_champaign_2025.csv", newline="", encoding="utf-8") as f:
        w = list(csv.DictReader(f))
    plant = date(2025, 4, 25)
    mid = [r for r in w if 71 <= (date.fromisoformat(r["date"]) - plant).days + 1 <= 120]
    u2 = sum(float(r["wind_2m_m_s"]) for r in mid) / len(mid)
    rh = sum(float(r["rh_min_pct"]) for r in mid) / len(mid)
    expected = 1.20 + (0.04 * (u2 - 2) - 0.004 * (rh - 45)) * (2.0 / 3.0) ** 0.3
    with open(tmp_path / "tables" / "daily_water_balance.csv", newline="", encoding="utf-8") as f:
        daily = list(csv.DictReader(f))
    assert float(daily[90]["kc"]) == pytest.approx(expected, abs=1e-3)
    rep = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert rep.index("## Warnings") < rep.index("## Set-up") and "tile-drained" in rep
