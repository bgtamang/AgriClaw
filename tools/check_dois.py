#!/usr/bin/env python3
"""Check every DOI cited in the skills against Crossref (or DataCite).

For each reference line that contains a DOI, the script looks the DOI up and
compares the registered first author and year with the line. A DOI that does
not resolve, or a line whose author or year disagrees with the registry, is
reported. Needs internet access; not part of CI.

Usage:
    python tools/check_dois.py                  # all skills
    python tools/check_dois.py claw-gdd ...     # selected skills
"""
from __future__ import annotations

import json
import re
import sys
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SKILLS = ROOT / "skills"
DOI_RE = re.compile(r"10\.\d{4,9}/[^\s<>\"'\]]+")
UA = "AgriClaw-check_dois/1.0 (https://github.com/bgtamang/AgriClaw)"


def clean_doi(raw: str) -> str:
    """Drop trailing punctuation and any unbalanced closing parenthesis.

    Elsevier DOIs such as 10.1016/0304-3800(94)00034-F contain parentheses,
    so they cannot simply be excluded from the match.
    """
    doi = raw.rstrip(".,;:")
    while doi.endswith(")") and doi.count(")") > doi.count("("):
        doi = doi[:-1].rstrip(".,;:")
    return doi


def fold(text: str) -> str:
    """Lower-case and strip accents so 'Salmerón' matches 'Salmeron'."""
    norm = unicodedata.normalize("NFKD", text)
    return "".join(c for c in norm if not unicodedata.combining(c)).lower()


def fetch_json(url: str, tries: int = 5) -> dict | None:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.load(resp)
        except urllib.error.HTTPError as err:
            if err.code == 404:
                return None
            if err.code in (429, 500, 502, 503, 504) and attempt < tries - 1:
                time.sleep(2 ** attempt)  # back off when rate-limited
                continue
            raise
    return None


def lookup(doi: str) -> dict:
    quoted = urllib.parse.quote(doi, safe="/")
    try:
        data = fetch_json(f"https://api.crossref.org/works/{quoted}")
        if data:
            msg = data["message"]
            authors = msg.get("author") or msg.get("editor") or []
            family = authors[0].get("family", authors[0].get("name", "")) if authors else ""
            parts = (msg.get("issued") or {}).get("date-parts") or [[None]]
            return {"ok": True, "registry": "Crossref",
                    "title": (msg.get("title") or [""])[0],
                    "family": family, "year": parts[0][0]}
        data = fetch_json(f"https://api.datacite.org/dois/{quoted}")
        if data:
            attrs = data["data"]["attributes"]
            creators = attrs.get("creators") or []
            family = ""
            if creators:
                family = creators[0].get("familyName") or creators[0].get("name", "").split(",")[0]
            # A concept DOI (e.g. Zenodo "all versions") takes the year of its
            # latest version, so its registered year is not checked.
            concept = any(r.get("relationType") == "HasVersion"
                          for r in attrs.get("relatedIdentifiers") or [])
            return {"ok": True, "registry": "DataCite",
                    "title": (attrs.get("titles") or [{}])[0].get("title", ""),
                    "family": family,
                    "year": None if concept else attrs.get("publicationYear")}
        return {"ok": False, "error": "not found in Crossref or DataCite"}
    except Exception as err:  # network trouble is reported, not fatal
        return {"ok": False, "error": f"lookup failed: {err}"}


def references(skill_dirs: list[Path]) -> dict[str, list[tuple[str, str]]]:
    """Map DOI -> [(skill, reference line), ...]."""
    found: dict[str, list[tuple[str, str]]] = {}
    for d in skill_dirs:
        md = d / "SKILL.md"
        if not md.exists():
            continue
        for line in md.read_text(encoding="utf-8").splitlines():
            for raw in DOI_RE.findall(line):
                found.setdefault(clean_doi(raw), []).append((d.name, line.strip()))
    return found


def main(argv: list[str]) -> int:
    names = argv[1:]
    dirs = [SKILLS / n for n in names] if names else sorted(p for p in SKILLS.iterdir() if p.is_dir())
    refs = references(dirs)
    if not refs:
        print("No DOIs found.")
        return 0
    results = {}
    for doi in refs:  # one at a time: Crossref rate-limits parallel clients
        results[doi] = lookup(doi)
        time.sleep(0.2)

    problems = 0
    for doi in sorted(refs):
        res = results[doi]
        for skill, line in refs[doi]:
            notes = []
            if not res["ok"]:
                notes.append(res["error"])
            else:
                fam = fold(res["family"] or "")
                if fam and fam not in fold(line):
                    notes.append(f"first author is {res['family']!r}")
                if res["year"] and str(res["year"]) not in line:
                    notes.append(f"registered year is {res['year']}")
            flag = "CHECK" if notes else "ok"
            problems += bool(notes)
            title = res.get("title", "")[:70] if res["ok"] else ""
            print(f"{flag:5}  {skill:24} {doi}  {title}" + (f"  <- {'; '.join(notes)}" if notes else ""))
    print(f"\n{len(refs)} DOIs, {problems} reference line(s) to check.")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
