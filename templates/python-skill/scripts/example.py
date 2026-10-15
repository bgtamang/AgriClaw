#!/usr/bin/env python3
"""claw-example: summarize a numeric trait by group (n, mean, SD, CV) and flag
groups with too few observations. Template for AgriClaw Python skills: copy
this file, keep the helper section unchanged, and replace the analysis section.

Usage:
    python scripts/example.py --demo
    python scripts/example.py --input data.csv --trait yield_kg_ha --group genotype

Part of AgriClaw <https://github.com/bgtamang/AgriClaw>. MIT license.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import math
import platform
import shlex
import statistics
import sys
from datetime import datetime
from pathlib import Path

SKILL = "claw-example"
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


# ---------------------------------------------------------------------------
# Analysis (replace in a real skill)
# ---------------------------------------------------------------------------


def read_rows(path: str) -> list[dict]:
    p = Path(path)
    if not p.is_file():
        raise SystemExit(f"Input file not found: {path}")
    with p.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def summarize_trait(rows: list[dict], trait: str, group: str, min_n: int = 3) -> list[dict]:
    if not rows:
        raise ValueError("Input has no rows.")
    missing = [c for c in (trait, group) if c not in rows[0]]
    if missing:
        raise ValueError(f"Column(s) not found: {', '.join(missing)}. "
                         f"Available: {', '.join(rows[0].keys())}")
    values: dict[str, list[float]] = {}
    for r in rows:
        raw = (r[trait] or "").strip()
        if raw in ("", "NA", "NaN"):
            continue
        try:
            x = float(raw)
        except ValueError:
            raise ValueError(f"Column '{trait}' must be numeric; found '{raw}'.") from None
        values.setdefault(r[group], []).append(x)
    out = []
    for g in sorted(values):
        x = values[g]
        n = len(x)
        mean = statistics.fmean(x)
        sd = statistics.stdev(x) if n > 1 else math.nan
        out.append({group: g, "n": n, "mean": mean, "sd": sd,
                    "cv_pct": 100 * sd / mean if mean else math.nan,
                    "low_n": n < min_n})
    return out


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(prog="example.py", description=__doc__.splitlines()[0])
    ap.add_argument("--demo", action="store_true", help="run on demo/example.csv")
    ap.add_argument("--input", help="CSV file, one row per observation")
    ap.add_argument("--trait", default="yield_kg_ha", help="numeric column (default: yield_kg_ha)")
    ap.add_argument("--group", default="genotype", help="grouping column (default: genotype)")
    ap.add_argument("--min-n", type=int, default=3, help="flag groups below this size (default: 3)")
    ap.add_argument("--out", default=f"agriclaw_output/{SKILL}", help="output folder")
    args = ap.parse_args(argv)

    if args.demo:
        args.input = str(SKILL_DIR / "demo" / "example.csv")
    if not args.input:
        ap.error("provide --input FILE.csv, or --demo")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    try:
        res = summarize_trait(read_rows(args.input), args.trait, args.group, args.min_n)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    save_table(res, out, "summary.csv")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.barh([r[args.group] for r in res], [r["mean"] for r in res],
            xerr=[0 if math.isnan(r["sd"]) else r["sd"] for r in res])
    ax.set_xlabel(args.trait)
    ax.set_title("Group means (± SD)")
    fig.tight_layout()
    (out / "figures").mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "figures" / "means.png", dpi=150)
    plt.close(fig)

    warn = [f"{sum(r['low_n'] for r in res)} group(s) have fewer than {args.min_n} observations."] \
        if any(r["low_n"] for r in res) else []
    lines = [f"# {SKILL} report", ""]
    if warn:
        lines += ["## Warnings", ""] + [f"- {w}" for w in warn] + [""]
    lines += ["## Summary", "", f"Input: `{args.input}`", "",
              "| group | n | mean | sd | cv (%) |", "|---|---|---|---|---|"]
    lines += [f"| {r[args.group]} | {r['n']} | {r['mean']:.4g} | {r['sd']:.3g} | {r['cv_pct']:.3g} |"
              for r in res]
    lines += ["", "![Group means](figures/means.png)", "", f"*{SKILL} {VERSION}*"]
    write_text(out / "report.md", "\n".join(lines) + "\n")
    write_repro_bundle(out, argv, inputs=[args.input], packages=["matplotlib"])
    print(f"Done: {len(res)} groups. Report: {(out / 'report.md').resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
