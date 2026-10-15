#!/usr/bin/env python3
"""Build catalog.json and the README skills table from every skills/*/SKILL.md.

Usage:
    python tools/build_catalog.py           # write catalog.json and update README.md
    python tools/build_catalog.py --check   # exit 1 if either file is out of date

The README table sits between the markers
<!-- skills:start --> and <!-- skills:end -->.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
FRONT_RE = re.compile(r"\A---\r?\n(.*?)\r?\n---", re.DOTALL)
DOMAIN_ORDER = ["meta", "field", "trials", "genomics", "transcriptomics", "remote-sensing", "physiology"]
DOMAIN_TITLES = {"meta": "Start here", "field": "Field and environment", "trials": "Field trials",
                 "genomics": "Genotypes and genetics", "transcriptomics": "Gene expression",
                 "remote-sensing": "Drone and satellite imagery", "physiology": "Plant physiology"}
STATUS_BADGE = {"tested": "tested", "experimental": "experimental", "guide": "guide"}


def load_skills() -> list[dict]:
    skills = []
    for skill_md in sorted((ROOT / "skills").glob("*/SKILL.md")):
        m = FRONT_RE.match(skill_md.read_text(encoding="utf-8"))
        if not m:
            continue
        fm = yaml.safe_load(m.group(1)) or {}
        meta = fm.get("metadata") or {}
        skills.append({
            "name": fm.get("name", skill_md.parent.name),
            "description": " ".join(str(fm.get("description", "")).split()),
            "summary": str(meta.get("summary", "")),
            "domain": meta.get("domain", "?"),
            "language": meta.get("language", "?"),
            "status": meta.get("status", "?"),
            "version": meta.get("version", "?"),
            "path": skill_md.parent.relative_to(ROOT).as_posix(),
        })
    order = {d: i for i, d in enumerate(DOMAIN_ORDER)}
    skills.sort(key=lambda s: (order.get(s["domain"], 99), s["name"]))
    return skills


def first_sentence(text: str, limit: int = 140) -> str:
    s = re.split(r"(?<=[.!?])\s", text, maxsplit=1)[0]
    return s if len(s) <= limit else s[: limit - 1].rstrip() + "…"


def readme_table(skills: list[dict]) -> str:
    lines = []
    for domain in DOMAIN_ORDER:
        group = [s for s in skills if s["domain"] == domain]
        if not group:
            continue
        lines += [f"### {DOMAIN_TITLES[domain]}", "", "| Skill | What it does | Language | Status |",
                  "|---|---|---|---|"]
        for s in group:
            what = s["summary"] or first_sentence(s["description"])
            lines.append(f"| [`{s['name']}`]({s['path']}/) | {what} "
                         f"| {s['language']} | {STATUS_BADGE.get(s['status'], s['status'])} |")
        lines.append("")
    counts = {k: sum(s["status"] == k for s in skills) for k in ("tested", "experimental", "guide")}
    lines.append(f"{len(skills)} skills: {counts['tested']} tested, {counts['experimental']} experimental, "
                 f"{counts['guide']} guide.")
    return "\n".join(lines)


def main(argv: list[str]) -> int:
    check = "--check" in argv
    skills = load_skills()
    catalog = json.dumps({"skills": skills}, indent=2, ensure_ascii=False) + "\n"
    cat_path = ROOT / "catalog.json"
    readme_path = ROOT / "README.md"
    readme = readme_path.read_text(encoding="utf-8") if readme_path.exists() else ""
    start, end = "<!-- skills:start -->", "<!-- skills:end -->"
    stale = []
    if start in readme and end in readme:
        new_readme = re.sub(re.escape(start) + r".*?" + re.escape(end),
                            lambda _: f"{start}\n{readme_table(skills)}\n{end}", readme, flags=re.DOTALL)
    else:
        new_readme = readme
        print("README.md has no skills markers; table not inserted.")
    if (cat_path.read_text(encoding="utf-8") if cat_path.exists() else "") != catalog:
        stale.append("catalog.json")
    if new_readme != readme:
        stale.append("README.md")
    if check:
        if stale:
            print("Out of date: " + ", ".join(stale) + ". Run python tools/build_catalog.py")
            return 1
        print("catalog.json and README.md are up to date.")
        return 0
    cat_path.write_text(catalog, encoding="utf-8")
    readme_path.write_text(new_readme, encoding="utf-8")
    print(f"Wrote catalog.json ({len(skills)} skills)" + ("; updated README.md" if "README.md" in stale else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
