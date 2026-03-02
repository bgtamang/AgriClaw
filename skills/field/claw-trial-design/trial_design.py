#!/usr/bin/env python3
"""
claw-trial-design: Field Trial Randomisation and Layout

Generates statistically valid experimental designs for crop breeding and
agronomy trials. Supports 9 design types with serpentine field layout,
publication-quality plot maps, and reproducibility bundles.

Supported designs:
    rcbd          Randomised Complete Block Design
    alpha         Alpha-lattice (incomplete block)
    augmented     Augmented (Federer) with replicated checks
    prep          Partially replicated
    split-plot    Split-plot (whole-plot + sub-plot)
    strip-plot    Strip-plot / split-block
    factorial     Factorial in RCBD
    crd           Completely Randomised Design
    latin-square  Latin Square

Usage:
    python trial_design.py --demo
    python trial_design.py --design rcbd --entries 24 --reps 4
    python trial_design.py --design alpha --entries 200 --reps 2
    python trial_design.py --design augmented --entries 500 --checks 5 --check-reps 3
    python trial_design.py --design split-plot --whole-plot 3 --sub-plot 4 --reps 4

References:
    RCBD, CRD, Latin Square:
        Cochran & Cox (1957). Experimental Designs, 2nd ed. Wiley.
    Alpha-lattice:
        Patterson & Williams (1976). A new class of resolvable incomplete
        block designs. Biometrika, 63(1), 83-92.
    Augmented design:
        Federer (1961). Augmented designs. Hawaiian Planters' Record, 56, 55-61.
    Partially replicated:
        Cullis et al. (2006). Analysis of yield trials using spatial methods.
        Journal of Agricultural Science, 144, 515-525.
    Split-plot:
        Steel & Torrie (1980). Principles and Procedures of Statistics,
        2nd ed. McGraw-Hill.
    Strip-plot:
        Gomez & Gomez (1984). Statistical Procedures for Agricultural
        Research, 2nd ed. Wiley.
    Factorial:
        Montgomery (2017). Design and Analysis of Experiments, 9th ed. Wiley.
"""

import argparse
import csv
import hashlib
import math
import os
import random
import sys
from datetime import datetime

# ---------------------------------------------------------------------------
# Shared utilities
# ---------------------------------------------------------------------------

def serpentine_layout(n_plots, n_rows, n_cols):
    """Generate serpentine (boustrophedon) plot numbering.

    Returns list of (row, col) tuples for plot indices 0..n_plots-1.
    Odd rows go left-to-right, even rows right-to-left (1-indexed rows).
    """
    positions = []
    plot = 0
    for row in range(1, n_rows + 1):
        if row % 2 == 1:
            cols = range(1, n_cols + 1)
        else:
            cols = range(n_cols, 0, -1)
        for col in cols:
            if plot < n_plots:
                positions.append((row, col))
                plot += 1
    return positions


def auto_field_dims(n_plots):
    """Choose reasonable row x col dimensions for a given number of plots."""
    sqrt_n = int(math.ceil(math.sqrt(n_plots)))
    # Try to find dimensions close to square
    for n_cols in range(sqrt_n, 0, -1):
        n_rows = math.ceil(n_plots / n_cols)
        if n_rows * n_cols >= n_plots:
            return n_rows, n_cols
    return n_plots, 1


def validate_design(field_book):
    """Basic validation of a field book."""
    issues = []
    # Check for missing plot numbers
    plots = [r["plot"] for r in field_book]
    if len(plots) != len(set(plots)):
        issues.append("Duplicate plot numbers detected")
    if not plots:
        issues.append("Empty field book")
    return issues


def write_csv(filepath, rows, fieldnames):
    """Write list of dicts to CSV."""
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


# ---------------------------------------------------------------------------
# Design generators
# Each returns (field_book, design_info)
#   field_book: list of dicts with at least {plot, row, range, block, entry, rep}
#   design_info: dict with metadata about the design
# ---------------------------------------------------------------------------

def design_rcbd(n_entries, n_reps, seed=None):
    """Randomised Complete Block Design.

    Reference: Cochran & Cox (1957). Experimental Designs, 2nd ed. Wiley.
    """
    rng = random.Random(seed)
    entries = list(range(1, n_entries + 1))
    n_plots = n_entries * n_reps
    n_rows, n_cols = auto_field_dims(n_plots)
    positions = serpentine_layout(n_plots, n_rows, n_cols)

    field_book = []
    plot_idx = 0
    for rep in range(1, n_reps + 1):
        block_entries = entries[:]
        rng.shuffle(block_entries)
        for entry in block_entries:
            row, col = positions[plot_idx]
            field_book.append({
                "plot": plot_idx + 1,
                "row": row,
                "range": col,
                "block": rep,
                "entry": f"E{entry:04d}",
                "rep": rep,
                "role": "test",
            })
            plot_idx += 1

    info = {
        "design": "RCBD",
        "entries": n_entries,
        "reps": n_reps,
        "total_plots": n_plots,
        "rows": n_rows,
        "cols": n_cols,
        "seed": seed,
    }
    return field_book, info


def design_alpha_lattice(n_entries, n_reps, block_size=None, seed=None):
    """Alpha-lattice (resolvable incomplete block design).

    Default block size: sqrt(n_entries), rounded to nearest valid divisor.

    Reference: Patterson & Williams (1976). A new class of resolvable
    incomplete block designs. Biometrika, 63(1), 83-92.
    """
    rng = random.Random(seed)

    # Find block size: nearest divisor of n_entries to sqrt(n_entries)
    if block_size is None:
        target = int(math.sqrt(n_entries))
        # Search for nearest divisor
        for offset in range(0, n_entries):
            if (target + offset) > 0 and n_entries % (target + offset) == 0:
                block_size = target + offset
                break
            if (target - offset) > 0 and n_entries % (target - offset) == 0:
                block_size = target - offset
                break
        if block_size is None or block_size < 2:
            block_size = n_entries  # Fallback to RCBD-like

    n_blocks_per_rep = n_entries // block_size
    n_plots = n_entries * n_reps
    n_rows, n_cols = auto_field_dims(n_plots)
    positions = serpentine_layout(n_plots, n_rows, n_cols)

    entries = list(range(1, n_entries + 1))
    field_book = []
    plot_idx = 0
    block_counter = 0

    for rep in range(1, n_reps + 1):
        rep_entries = entries[:]
        rng.shuffle(rep_entries)

        for b in range(n_blocks_per_rep):
            block_counter += 1
            block_entries = rep_entries[b * block_size: (b + 1) * block_size]
            rng.shuffle(block_entries)
            for entry in block_entries:
                row, col = positions[plot_idx]
                field_book.append({
                    "plot": plot_idx + 1,
                    "row": row,
                    "range": col,
                    "block": block_counter,
                    "entry": f"E{entry:04d}",
                    "rep": rep,
                    "role": "test",
                })
                plot_idx += 1

    info = {
        "design": "Alpha-Lattice",
        "entries": n_entries,
        "reps": n_reps,
        "block_size": block_size,
        "blocks_per_rep": n_blocks_per_rep,
        "total_blocks": block_counter,
        "total_plots": n_plots,
        "rows": n_rows,
        "cols": n_cols,
        "seed": seed,
    }
    return field_book, info


def design_augmented(n_test, n_checks, check_reps=3, seed=None):
    """Augmented (Federer) design.

    Test entries unreplicated, checks replicated systematically.
    Checks placed every k-th plot for even spatial coverage.

    Reference: Federer (1961). Augmented designs.
    Hawaiian Planters' Record, 56, 55-61.
    """
    rng = random.Random(seed)

    total_check_plots = n_checks * check_reps
    n_plots = n_test + total_check_plots

    # Check spacing: place a check every k-th plot
    k = max(1, n_plots // total_check_plots)

    n_rows, n_cols = auto_field_dims(n_plots)
    positions = serpentine_layout(n_plots, n_rows, n_cols)

    # Build plot assignment
    test_entries = list(range(1, n_test + 1))
    rng.shuffle(test_entries)

    # Create check pool: each check repeated check_reps times
    check_pool = []
    for rep in range(check_reps):
        for c in range(1, n_checks + 1):
            check_pool.append((c, rep + 1))
    rng.shuffle(check_pool)

    # Assign plots: systematic check placement
    field_book = []
    test_idx = 0
    check_idx = 0
    block = 1

    for plot_num in range(n_plots):
        row, col = positions[plot_num]

        # Place check at every k-th plot (and ensure we use all checks)
        if (plot_num % k == 0 or n_plots - plot_num <= len(check_pool) - check_idx) and check_idx < len(check_pool):
            c_id, c_rep = check_pool[check_idx]
            field_book.append({
                "plot": plot_num + 1,
                "row": row,
                "range": col,
                "block": block,
                "entry": f"CHK{c_id:02d}",
                "rep": c_rep,
                "role": "check",
            })
            check_idx += 1
        else:
            if test_idx < len(test_entries):
                field_book.append({
                    "plot": plot_num + 1,
                    "row": row,
                    "range": col,
                    "block": block,
                    "entry": f"E{test_entries[test_idx]:04d}",
                    "rep": 1,
                    "role": "test",
                })
                test_idx += 1

        # Block boundary: every n_cols plots (one row = one block)
        if (plot_num + 1) % n_cols == 0:
            block += 1

    info = {
        "design": "Augmented (Federer)",
        "test_entries": n_test,
        "checks": n_checks,
        "check_reps": check_reps,
        "total_plots": n_plots,
        "check_spacing": k,
        "rows": n_rows,
        "cols": n_cols,
        "seed": seed,
    }
    return field_book, info


def design_prep(n_entries, n_reps=2, rep_fraction=0.3, seed=None):
    """Partially replicated (p-rep) design.

    A fraction of entries get full replication, the rest get 1 rep.

    Reference: Cullis et al. (2006). Analysis of yield trials using
    spatial methods. J. Agric. Sci., 144, 515-525.
    """
    rng = random.Random(seed)

    n_replicated = max(1, int(n_entries * rep_fraction))
    n_unreplicated = n_entries - n_replicated
    n_plots = n_replicated * n_reps + n_unreplicated

    n_rows, n_cols = auto_field_dims(n_plots)
    positions = serpentine_layout(n_plots, n_rows, n_cols)

    entries = list(range(1, n_entries + 1))
    rng.shuffle(entries)

    replicated = entries[:n_replicated]
    unreplicated = entries[n_replicated:]

    # Build all plots
    all_plots = []
    for entry in replicated:
        for rep in range(1, n_reps + 1):
            all_plots.append((entry, rep, "replicated"))
    for entry in unreplicated:
        all_plots.append((entry, 1, "unreplicated"))

    rng.shuffle(all_plots)

    field_book = []
    block = 1
    for i, (entry, rep, role) in enumerate(all_plots):
        row, col = positions[i]
        field_book.append({
            "plot": i + 1,
            "row": row,
            "range": col,
            "block": block,
            "entry": f"E{entry:04d}",
            "rep": rep,
            "role": role,
        })
        if (i + 1) % n_cols == 0:
            block += 1

    info = {
        "design": "Partially Replicated (p-rep)",
        "entries": n_entries,
        "replicated_entries": n_replicated,
        "unreplicated_entries": n_unreplicated,
        "reps_for_replicated": n_reps,
        "rep_fraction": rep_fraction,
        "total_plots": n_plots,
        "rows": n_rows,
        "cols": n_cols,
        "seed": seed,
    }
    return field_book, info


def design_split_plot(n_whole, n_sub, n_reps, seed=None):
    """Split-plot design.

    Whole-plot factor randomised to main plots within each rep.
    Sub-plot factor randomised within each whole-plot.

    Expert decision: whole-plot randomised FIRST, then sub-plot within.
    Common error: randomising both independently.

    Reference: Steel & Torrie (1980). Principles and Procedures of
    Statistics, 2nd ed. McGraw-Hill.
    """
    rng = random.Random(seed)

    n_plots = n_whole * n_sub * n_reps
    n_rows, n_cols = auto_field_dims(n_plots)
    positions = serpentine_layout(n_plots, n_rows, n_cols)

    field_book = []
    plot_idx = 0
    block = 0

    for rep in range(1, n_reps + 1):
        whole_order = list(range(1, n_whole + 1))
        rng.shuffle(whole_order)  # Randomise whole-plot factor

        for wp in whole_order:
            block += 1
            sub_order = list(range(1, n_sub + 1))
            rng.shuffle(sub_order)  # Randomise sub-plot within whole-plot

            for sp in sub_order:
                row, col = positions[plot_idx]
                field_book.append({
                    "plot": plot_idx + 1,
                    "row": row,
                    "range": col,
                    "block": block,
                    "entry": f"W{wp}_S{sp}",
                    "rep": rep,
                    "role": f"whole={wp} sub={sp}",
                    "whole_plot": wp,
                    "sub_plot": sp,
                })
                plot_idx += 1

    info = {
        "design": "Split-Plot",
        "whole_plot_levels": n_whole,
        "sub_plot_levels": n_sub,
        "reps": n_reps,
        "total_plots": n_plots,
        "rows": n_rows,
        "cols": n_cols,
        "seed": seed,
    }
    return field_book, info


def design_strip_plot(n_row_factor, n_col_factor, n_reps, seed=None):
    """Strip-plot (split-block) design.

    Row factor applied in horizontal strips, column factor in vertical strips.
    Each rep is a row_factor x col_factor grid.

    Reference: Gomez & Gomez (1984). Statistical Procedures for
    Agricultural Research, 2nd ed. Wiley.
    """
    rng = random.Random(seed)

    plots_per_rep = n_row_factor * n_col_factor
    n_plots = plots_per_rep * n_reps

    # Layout: each rep is n_row_factor rows x n_col_factor cols
    total_rows = n_row_factor * n_reps
    total_cols = n_col_factor

    positions = serpentine_layout(n_plots, total_rows, total_cols)

    field_book = []
    plot_idx = 0

    for rep in range(1, n_reps + 1):
        row_order = list(range(1, n_row_factor + 1))
        col_order = list(range(1, n_col_factor + 1))
        rng.shuffle(row_order)
        rng.shuffle(col_order)

        for ri, rf in enumerate(row_order):
            for ci, cf in enumerate(col_order):
                row, col = positions[plot_idx]
                field_book.append({
                    "plot": plot_idx + 1,
                    "row": row,
                    "range": col,
                    "block": rep,
                    "entry": f"R{rf}_C{cf}",
                    "rep": rep,
                    "role": f"row_factor={rf} col_factor={cf}",
                    "row_factor": rf,
                    "col_factor": cf,
                })
                plot_idx += 1

    info = {
        "design": "Strip-Plot (Split-Block)",
        "row_factor_levels": n_row_factor,
        "col_factor_levels": n_col_factor,
        "reps": n_reps,
        "total_plots": n_plots,
        "rows": total_rows,
        "cols": total_cols,
        "seed": seed,
    }
    return field_book, info


def design_factorial(factor_levels, n_reps, seed=None):
    """Factorial design in RCBD.

    All combinations of factor levels, arranged in randomised blocks.

    Reference: Montgomery (2017). Design and Analysis of Experiments,
    9th ed. Wiley.
    """
    rng = random.Random(seed)

    # Generate all treatment combinations
    def cartesian(lists):
        if not lists:
            return [()]
        result = []
        for item in lists[0]:
            for rest in cartesian(lists[1:]):
                result.append((item,) + rest)
        return result

    factor_ranges = [list(range(1, fl + 1)) for fl in factor_levels]
    treatments = cartesian(factor_ranges)
    n_treatments = len(treatments)
    n_plots = n_treatments * n_reps

    n_rows, n_cols = auto_field_dims(n_plots)
    positions = serpentine_layout(n_plots, n_rows, n_cols)

    field_book = []
    plot_idx = 0

    for rep in range(1, n_reps + 1):
        rep_treatments = list(range(n_treatments))
        rng.shuffle(rep_treatments)

        for ti in rep_treatments:
            trt = treatments[ti]
            row, col = positions[plot_idx]
            label = "x".join(f"F{i+1}L{v}" for i, v in enumerate(trt))
            field_book.append({
                "plot": plot_idx + 1,
                "row": row,
                "range": col,
                "block": rep,
                "entry": label,
                "rep": rep,
                "role": "treatment",
            })
            plot_idx += 1

    info = {
        "design": "Factorial in RCBD",
        "factor_levels": factor_levels,
        "n_factors": len(factor_levels),
        "treatments": n_treatments,
        "reps": n_reps,
        "total_plots": n_plots,
        "rows": n_rows,
        "cols": n_cols,
        "seed": seed,
    }
    return field_book, info


def design_crd(n_entries, n_reps, seed=None):
    """Completely Randomised Design.

    No blocking. All plots fully randomised.

    Reference: Cochran & Cox (1957). Experimental Designs, 2nd ed. Wiley.
    """
    rng = random.Random(seed)

    n_plots = n_entries * n_reps
    n_rows, n_cols = auto_field_dims(n_plots)
    positions = serpentine_layout(n_plots, n_rows, n_cols)

    all_plots = []
    for entry in range(1, n_entries + 1):
        for rep in range(1, n_reps + 1):
            all_plots.append((entry, rep))

    rng.shuffle(all_plots)

    field_book = []
    for i, (entry, rep) in enumerate(all_plots):
        row, col = positions[i]
        field_book.append({
            "plot": i + 1,
            "row": row,
            "range": col,
            "block": 1,
            "entry": f"E{entry:04d}",
            "rep": rep,
            "role": "test",
        })

    info = {
        "design": "CRD",
        "entries": n_entries,
        "reps": n_reps,
        "total_plots": n_plots,
        "rows": n_rows,
        "cols": n_cols,
        "seed": seed,
    }
    return field_book, info


def design_latin_square(n_entries, seed=None):
    """Latin Square design.

    n x n grid where each entry appears exactly once per row and column.
    Two blocking factors (row and column).

    Reference: Cochran & Cox (1957). Experimental Designs, 2nd ed. Wiley.
    """
    if n_entries > 12:
        print(f"Warning: Latin square with n={n_entries} is very large. "
              f"Consider n <= 12.", file=sys.stderr)

    rng = random.Random(seed)
    n = n_entries

    # Generate a standard Latin square by cyclic permutation, then shuffle
    square = []
    for i in range(n):
        row = [(i + j) % n + 1 for j in range(n)]
        square.append(row)

    # Shuffle rows
    rng.shuffle(square)
    # Shuffle columns
    cols_order = list(range(n))
    rng.shuffle(cols_order)
    square = [[row[c] for c in cols_order] for row in square]

    field_book = []
    plot = 0
    for r in range(n):
        for c in range(n):
            plot += 1
            field_book.append({
                "plot": plot,
                "row": r + 1,
                "range": c + 1,
                "block": r + 1,
                "entry": f"E{square[r][c]:04d}",
                "rep": 1,
                "role": "test",
            })

    info = {
        "design": "Latin Square",
        "entries": n_entries,
        "size": f"{n}x{n}",
        "total_plots": n * n,
        "rows": n,
        "cols": n,
        "seed": seed,
    }
    return field_book, info


# ---------------------------------------------------------------------------
# Output generation
# ---------------------------------------------------------------------------

def generate_field_layout(field_book, design_info, output_dir):
    """Generate field layout plot with color-coded entries."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
    except ImportError:
        print("matplotlib not installed -- skipping plot.", file=sys.stderr)
        return False

    fig_dir = os.path.join(output_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    n_rows = design_info["rows"]
    n_cols = design_info["cols"]

    # Build grid
    grid = {}
    entries_set = set()
    for p in field_book:
        grid[(p["row"], p["range"])] = p
        entries_set.add(p["entry"])

    # Color mapping
    entries_list = sorted(entries_set)
    n_entries = len(entries_list)

    # Use checks as red, test entries as color spectrum
    check_entries = {p["entry"] for p in field_book if p["role"] == "check"}

    if n_entries <= 20:
        cmap = plt.cm.tab20
    else:
        cmap = plt.cm.nipy_spectral

    entry_colors = {}
    test_entries = [e for e in entries_list if e not in check_entries]
    for i, e in enumerate(test_entries):
        entry_colors[e] = cmap(i / max(1, len(test_entries) - 1))
    for e in check_entries:
        entry_colors[e] = (1.0, 0.2, 0.2, 1.0)  # Red for checks

    # Figure size scales with field dimensions
    fig_width = max(8, n_cols * 0.6)
    fig_height = max(4, n_rows * 0.5)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    for r in range(1, n_rows + 1):
        for c in range(1, n_cols + 1):
            p = grid.get((r, c))
            if p:
                color = entry_colors.get(p["entry"], (0.8, 0.8, 0.8, 1.0))
                rect = plt.Rectangle((c - 0.5, r - 0.5), 1, 1,
                                     facecolor=color, edgecolor="white",
                                     linewidth=0.5)
                ax.add_patch(rect)

                # Label: entry name (truncated if too long)
                label = p["entry"]
                if len(label) > 8:
                    label = label[:7] + ".."
                fontsize = max(4, min(8, 120 // max(n_rows, n_cols)))
                ax.text(c, r, label, ha="center", va="center",
                        fontsize=fontsize, color="black")

    ax.set_xlim(0.5, n_cols + 0.5)
    ax.set_ylim(n_rows + 0.5, 0.5)
    ax.set_xlabel("Range")
    ax.set_ylabel("Row")
    ax.set_title(f"Field Layout: {design_info['design']} ({design_info['total_plots']} plots)")
    ax.set_aspect("equal")

    # Grid ticks
    ax.set_xticks(range(1, n_cols + 1))
    ax.set_yticks(range(1, n_rows + 1))

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "field_layout.png"), dpi=150)
    plt.close()
    return True


def generate_report(field_book, design_info, output_dir):
    """Generate markdown report."""
    n_checks = sum(1 for p in field_book if p["role"] == "check")
    n_test = sum(1 for p in field_book if p["role"] != "check")

    lines = [
        f"# Trial Design Report: {design_info['design']}",
        "",
        "## Design Parameters",
        "",
        "| Parameter | Value |",
        "|-----------|-------|",
    ]

    for key, val in design_info.items():
        if key == "seed" and val is None:
            val = "random"
        lines.append(f"| {key.replace('_', ' ').title()} | {val} |")

    lines += [
        "",
        "## Plot Summary",
        "",
        f"| Metric | Count |",
        f"|--------|-------|",
        f"| Total plots | {len(field_book)} |",
        f"| Field dimensions | {design_info['rows']} rows x {design_info['cols']} ranges |",
        f"| Test plots | {n_test} |",
        f"| Check plots | {n_checks} |",
        "",
        "## Field Layout",
        "",
        "![Field Layout](figures/field_layout.png)",
        "",
        "## Files",
        "",
        "- `tables/field_book.csv` -- Full field book ready for data collection",
        "- `tables/design_summary.csv` -- Design parameters",
        "- `figures/field_layout.png` -- Visual field map",
        "",
        "## References",
        "",
    ]

    # Design-specific references
    refs = {
        "RCBD": "Cochran & Cox (1957). Experimental Designs, 2nd ed. Wiley.",
        "Alpha-Lattice": "Patterson & Williams (1976). A new class of resolvable incomplete block designs. Biometrika, 63(1), 83-92.",
        "Augmented (Federer)": "Federer (1961). Augmented designs. Hawaiian Planters' Record, 56, 55-61.",
        "Partially Replicated (p-rep)": "Cullis et al. (2006). Analysis of yield trials using spatial methods. J. Agric. Sci., 144, 515-525.",
        "Split-Plot": "Steel & Torrie (1980). Principles and Procedures of Statistics, 2nd ed. McGraw-Hill.",
        "Strip-Plot (Split-Block)": "Gomez & Gomez (1984). Statistical Procedures for Agricultural Research, 2nd ed. Wiley.",
        "Factorial in RCBD": "Montgomery (2017). Design and Analysis of Experiments, 9th ed. Wiley.",
        "CRD": "Cochran & Cox (1957). Experimental Designs, 2nd ed. Wiley.",
        "Latin Square": "Cochran & Cox (1957). Experimental Designs, 2nd ed. Wiley.",
    }
    ref = refs.get(design_info["design"], "")
    if ref:
        lines.append(f"- {ref}")
    lines += [
        "",
        "---",
        f"*Generated by [AgriClaw claw-trial-design](https://github.com/bgtamang/AgriClaw) v0.1.0*",
    ]

    report_path = os.path.join(output_dir, "report.md")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def generate_repro_bundle(output_dir, cmd_str):
    """Generate reproducibility bundle."""
    repro_dir = os.path.join(output_dir, "reproducibility")
    os.makedirs(repro_dir, exist_ok=True)

    with open(os.path.join(repro_dir, "commands.sh"), "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"# Reproduce this analysis\n")
        f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"pip install pandas matplotlib\n")
        f.write(f"{cmd_str}\n")

    with open(os.path.join(repro_dir, "environment.yml"), "w") as f:
        f.write("name: claw-trial-design\n")
        f.write("dependencies:\n")
        f.write("  - python>=3.9\n")
        f.write("  - pip:\n")
        f.write("    - pandas\n")
        f.write("    - matplotlib\n")

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

DESIGNS = ["rcbd", "alpha", "augmented", "prep", "split-plot",
           "strip-plot", "factorial", "crd", "latin-square"]


def main():
    parser = argparse.ArgumentParser(
        description="claw-trial-design: Field Trial Randomisation (AgriClaw)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  python trial_design.py --demo\n"
               "  python trial_design.py --design rcbd --entries 24 --reps 4\n"
               "  python trial_design.py --design alpha --entries 200 --reps 2\n"
               "  python trial_design.py --design augmented --entries 500 --checks 5\n"
               "  python trial_design.py --design split-plot --whole-plot 3 --sub-plot 4 --reps 4\n"
               "  python trial_design.py --design factorial --factors 3,4,2 --reps 3\n"
               "  python trial_design.py --design latin-square --entries 5\n",
    )
    parser.add_argument("--demo", action="store_true",
                        help="Run demo: augmented design, 200 entries, 5 checks")
    parser.add_argument("--design", choices=DESIGNS,
                        help="Design type")
    parser.add_argument("--entries", type=int, help="Number of entries/treatments")
    parser.add_argument("--reps", type=int, help="Number of replications")
    parser.add_argument("--checks", type=int, default=5,
                        help="Number of check entries (augmented design)")
    parser.add_argument("--check-reps", type=int, default=3,
                        help="Reps per check (augmented design)")
    parser.add_argument("--block-size", type=int, default=None,
                        help="Block size (alpha-lattice; default: auto)")
    parser.add_argument("--rep-fraction", type=float, default=0.3,
                        help="Fraction of entries replicated (p-rep; default: 0.3)")
    parser.add_argument("--whole-plot", type=int,
                        help="Whole-plot factor levels (split-plot)")
    parser.add_argument("--sub-plot", type=int,
                        help="Sub-plot factor levels (split-plot)")
    parser.add_argument("--row-factor", type=int,
                        help="Row factor levels (strip-plot)")
    parser.add_argument("--col-factor", type=int,
                        help="Column factor levels (strip-plot)")
    parser.add_argument("--factors", type=str,
                        help="Factor levels comma-separated (factorial; e.g. 3,4,2)")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed for reproducibility")
    parser.add_argument("--output", default="output",
                        help="Output directory (default: output)")
    args = parser.parse_args()

    # Demo mode
    if args.demo:
        args.design = "augmented"
        args.entries = 200
        args.checks = 5
        args.check_reps = 3
        args.output = "demo_output"
        args.seed = 42
        print("Running demo: Augmented design, 200 test entries, 5 checks, 3 check reps")

    if not args.demo and not args.design:
        parser.error("Provide --design or --demo")

    # Generate seed if not provided
    if args.seed is None:
        args.seed = random.randint(1, 999999)
    print(f"Seed: {args.seed}")

    # Build command string for reproducibility
    cmd_parts = ["python trial_design.py", f"--design {args.design}", f"--seed {args.seed}"]

    # Dispatch to design generator
    if args.design == "rcbd":
        if not args.entries or not args.reps:
            parser.error("RCBD requires --entries and --reps")
        cmd_parts += [f"--entries {args.entries}", f"--reps {args.reps}"]
        field_book, info = design_rcbd(args.entries, args.reps, args.seed)

    elif args.design == "alpha":
        if not args.entries or not args.reps:
            parser.error("Alpha-lattice requires --entries and --reps")
        cmd_parts += [f"--entries {args.entries}", f"--reps {args.reps}"]
        if args.block_size:
            cmd_parts.append(f"--block-size {args.block_size}")
        field_book, info = design_alpha_lattice(args.entries, args.reps, args.block_size, args.seed)

    elif args.design == "augmented":
        if not args.entries:
            parser.error("Augmented requires --entries")
        cmd_parts += [f"--entries {args.entries}", f"--checks {args.checks}",
                      f"--check-reps {args.check_reps}"]
        field_book, info = design_augmented(args.entries, args.checks, args.check_reps, args.seed)

    elif args.design == "prep":
        if not args.entries:
            parser.error("P-rep requires --entries")
        reps = args.reps or 2
        cmd_parts += [f"--entries {args.entries}", f"--reps {reps}",
                      f"--rep-fraction {args.rep_fraction}"]
        field_book, info = design_prep(args.entries, reps, args.rep_fraction, args.seed)

    elif args.design == "split-plot":
        if not args.whole_plot or not args.sub_plot or not args.reps:
            parser.error("Split-plot requires --whole-plot, --sub-plot, and --reps")
        cmd_parts += [f"--whole-plot {args.whole_plot}", f"--sub-plot {args.sub_plot}",
                      f"--reps {args.reps}"]
        field_book, info = design_split_plot(args.whole_plot, args.sub_plot, args.reps, args.seed)

    elif args.design == "strip-plot":
        if not args.row_factor or not args.col_factor or not args.reps:
            parser.error("Strip-plot requires --row-factor, --col-factor, and --reps")
        cmd_parts += [f"--row-factor {args.row_factor}", f"--col-factor {args.col_factor}",
                      f"--reps {args.reps}"]
        field_book, info = design_strip_plot(args.row_factor, args.col_factor, args.reps, args.seed)

    elif args.design == "factorial":
        if not args.factors or not args.reps:
            parser.error("Factorial requires --factors and --reps")
        factor_levels = [int(x) for x in args.factors.split(",")]
        cmd_parts += [f"--factors {args.factors}", f"--reps {args.reps}"]
        field_book, info = design_factorial(factor_levels, args.reps, args.seed)

    elif args.design == "crd":
        if not args.entries or not args.reps:
            parser.error("CRD requires --entries and --reps")
        cmd_parts += [f"--entries {args.entries}", f"--reps {args.reps}"]
        field_book, info = design_crd(args.entries, args.reps, args.seed)

    elif args.design == "latin-square":
        if not args.entries:
            parser.error("Latin square requires --entries")
        cmd_parts += [f"--entries {args.entries}"]
        field_book, info = design_latin_square(args.entries, args.seed)

    # Validate
    issues = validate_design(field_book)
    if issues:
        for issue in issues:
            print(f"VALIDATION ERROR: {issue}", file=sys.stderr)
        sys.exit(1)

    # Output
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    tables_dir = os.path.join(output_dir, "tables")
    os.makedirs(tables_dir, exist_ok=True)

    # Determine field book columns based on design
    base_fields = ["plot", "row", "range", "block", "entry", "rep", "role"]
    extra_fields = []
    if args.design == "split-plot":
        extra_fields = ["whole_plot", "sub_plot"]
    elif args.design == "strip-plot":
        extra_fields = ["row_factor", "col_factor"]

    write_csv(os.path.join(tables_dir, "field_book.csv"), field_book, base_fields + extra_fields)
    print(f"Wrote field_book.csv ({len(field_book)} plots)")

    # Design summary
    summary_rows = [{"parameter": k, "value": str(v)} for k, v in info.items()]
    write_csv(os.path.join(tables_dir, "design_summary.csv"), summary_rows, ["parameter", "value"])

    # Plot
    plot_ok = generate_field_layout(field_book, info, output_dir)
    if plot_ok:
        print("Generated field_layout.png")

    # Report
    generate_report(field_book, info, output_dir)
    print("Generated report.md")

    # Reproducibility
    cmd_str = " ".join(cmd_parts)
    generate_repro_bundle(output_dir, cmd_str)
    print("Generated reproducibility bundle")

    # Summary
    print(f"\n{'='*50}")
    print(f"Design: {info['design']}")
    print(f"Total plots: {info['total_plots']}")
    print(f"Field: {info['rows']} rows x {info['cols']} ranges")
    print(f"Seed: {args.seed}")
    print(f"\nOutput: {os.path.abspath(output_dir)}/")


if __name__ == "__main__":
    main()
