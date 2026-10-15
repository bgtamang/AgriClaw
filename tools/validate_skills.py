#!/usr/bin/env python3
"""Check every skill folder against the Agent Skills specification
(https://agentskills.io/specification) and the AgriClaw conventions in
CONVENTIONS.md.

Usage:
    python tools/validate_skills.py                 # all skills/*
    python tools/validate_skills.py claw-gdd ...    # selected skills
    python tools/validate_skills.py --templates     # also check templates/*

Exits with status 1 if any skill has an error. Warnings don't fail the run.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]

ALLOWED_KEYS = {"name", "description", "license", "compatibility", "metadata", "allowed-tools"}
REQUIRED_META = ("version", "domain", "language", "status", "summary")
DOMAINS = {"field", "trials", "genomics", "transcriptomics", "remote-sensing", "physiology", "meta"}
STATUSES = {"tested", "experimental", "guide"}
LANGUAGES = {"R", "Python", "R+Python", "none"}
NAME_RE = re.compile(r"^[a-z0-9]+(-[a-z0-9]+)*$")
VERSION_RE = re.compile(r"^\d+\.\d+\.\d+$")
FRONT_RE = re.compile(r"\A---\r?\n(.*?)\r?\n---[ \t]*\r?\n(.*)\Z", re.DOTALL)

# Body sections every runnable skill needs (matched case-insensitively in headings).
SECTIONS_RUNNABLE = ["what it does", "when to use", "inputs", "how to run", "outputs",
                     "expert decisions", "checks", "references"]
SECTIONS_GUIDE = ["what it does", "when to use", "references"]


def check_skill(skill_dir: Path, is_template: bool = False) -> tuple[list[str], list[str], dict]:
    errors: list[str] = []
    warnings: list[str] = []
    info: dict = {"name": skill_dir.name, "status": "?", "domain": "?", "language": "?"}

    skill_md = skill_dir / "SKILL.md"
    if not skill_md.is_file():
        return ["SKILL.md is missing"], warnings, info
    text = skill_md.read_text(encoding="utf-8")
    m = FRONT_RE.match(text)
    if not m:
        return ["SKILL.md must start with a '---' YAML frontmatter block closed by '---'"], warnings, info
    try:
        fm = yaml.safe_load(m.group(1))
    except yaml.YAMLError as e:
        return [f"frontmatter is not valid YAML: {e}"], warnings, info
    if not isinstance(fm, dict):
        return ["frontmatter must be a YAML mapping"], warnings, info
    body = m.group(2)

    # --- top-level keys ---------------------------------------------------
    extra = sorted(set(fm) - ALLOWED_KEYS)
    if extra:
        errors.append(f"keys not allowed by the Agent Skills spec: {', '.join(extra)} "
                      f"(move them under metadata)")

    name = fm.get("name")
    if not isinstance(name, str) or not name:
        errors.append("name is required")
    else:
        if len(name) > 64:
            errors.append(f"name is {len(name)} characters (max 64)")
        if not NAME_RE.match(name):
            errors.append("name must be lowercase letters, digits and single hyphens")
        if not is_template and name != skill_dir.name:
            errors.append(f"name '{name}' must equal the folder name '{skill_dir.name}'")
        if not name.startswith("claw-") and name != "agriclaw":
            warnings.append("AgriClaw skill names start with 'claw-'")
        info["name"] = name

    desc = fm.get("description")
    if not isinstance(desc, str) or not desc.strip():
        errors.append("description is required")
    else:
        if len(desc) > 1024:
            errors.append(f"description is {len(desc)} characters (max 1024)")
        if "use when" not in desc.lower() and "use for" not in desc.lower() and "use as" not in desc.lower():
            warnings.append("description should say when to use the skill ('Use when ...')")

    if fm.get("license") != "MIT":
        warnings.append("license should be MIT")

    comp = fm.get("compatibility")
    if not isinstance(comp, str) or not comp.strip():
        errors.append("compatibility is required in AgriClaw (language, packages, internet needs)")
    elif len(comp) > 500:
        errors.append(f"compatibility is {len(comp)} characters (max 500)")

    # --- metadata -------------------------------------------------------------
    meta = fm.get("metadata")
    if not isinstance(meta, dict):
        errors.append("metadata must be a mapping with version, domain, language, status, summary")
        meta = {}
    for k, v in meta.items():
        if not isinstance(k, str) or not isinstance(v, str):
            errors.append(f"metadata.{k} must be a string (quote it in YAML); got {type(v).__name__}")
    for k in REQUIRED_META:
        if k not in meta:
            errors.append(f"metadata.{k} is required")
    status = str(meta.get("status", "?"))
    info.update(status=status, domain=str(meta.get("domain", "?")),
                language=str(meta.get("language", "?")), version=str(meta.get("version", "?")))
    if "version" in meta and not VERSION_RE.match(str(meta["version"])):
        errors.append(f"metadata.version '{meta['version']}' must look like 1.0.0")
    if "domain" in meta and meta["domain"] not in DOMAINS:
        errors.append(f"metadata.domain must be one of {sorted(DOMAINS)}")
    if "status" in meta and status not in STATUSES:
        errors.append(f"metadata.status must be one of {sorted(STATUSES)}")
    if "language" in meta and meta["language"] not in LANGUAGES:
        errors.append(f"metadata.language must be one of {sorted(LANGUAGES)}")
    if "summary" in meta and len(str(meta["summary"])) > 100:
        errors.append(f"metadata.summary is {len(str(meta['summary']))} characters (max 100)")

    # --- body -----------------------------------------------------------------
    n_lines = body.count("\n") + 1
    if n_lines > 500:
        errors.append(f"SKILL.md body is {n_lines} lines (spec maximum 500); move detail to references/")
    elif n_lines > 300:
        warnings.append(f"SKILL.md body is {n_lines} lines; AgriClaw aims for under 300")
    headings = [h.lower() for h in re.findall(r"^#{1,4}\s+(.+?)\s*$", body, flags=re.MULTILINE)]
    needed = SECTIONS_GUIDE if status == "guide" else SECTIONS_RUNNABLE
    for sec in needed:
        if not any(sec in h for h in headings):
            errors.append(f"SKILL.md has no section heading containing '{sec}'")
    ref_text = body.lower().split("references", 1)[-1] if "references" in body.lower() else ""
    if status != "guide" and not re.search(r"doi\.org/|doi:|https?://", ref_text):
        warnings.append("References section has no DOI or URL")

    # --- files ------------------------------------------------------------------
    scripts = [p for p in (skill_dir / "scripts").glob("*") if p.suffix in (".R", ".py")] \
        if (skill_dir / "scripts").is_dir() else []
    tests = [p for p in (skill_dir / "tests").glob("test_*") if p.suffix in (".R", ".py")] \
        if (skill_dir / "tests").is_dir() else []
    demo = [p for p in (skill_dir / "demo").rglob("*") if p.is_file()] \
        if (skill_dir / "demo").is_dir() else []
    if status in ("tested", "experimental") and not scripts:
        errors.append("scripts/ needs at least one .R or .py file")
    if status == "tested":
        if not tests:
            errors.append("status 'tested' needs tests/test_*.R or tests/test_*.py")
        if not demo:
            warnings.append("demo/ is empty: fine only if the demo data come from an R package via system.file()")
    for p in skill_dir.rglob("*"):
        if p.is_file() and p.stat().st_size > 1_000_000:
            warnings.append(f"{p.relative_to(skill_dir).as_posix()} is {p.stat().st_size / 1e6:.1f} MB (keep demo data small)")
    return errors, warnings, info


def main(argv: list[str]) -> int:
    with_templates = "--templates" in argv
    names = [a for a in argv if not a.startswith("--")]
    dirs = sorted(p for p in (ROOT / "skills").iterdir() if p.is_dir()) if (ROOT / "skills").is_dir() else []
    if names:
        dirs = [d for d in dirs if d.name in names]
        missing = set(names) - {d.name for d in dirs}
        if missing:
            print(f"No such skill(s): {', '.join(sorted(missing))}")
            return 1
    targets = [(d, False) for d in dirs]
    if with_templates:
        targets += [(d, True) for d in sorted((ROOT / "templates").iterdir()) if d.is_dir()]

    n_err = 0
    print(f"{'skill':28} {'status':13} {'domain':16} {'lang':9} result")
    for d, is_template in targets:
        errors, warnings, info = check_skill(d, is_template)
        n_err += bool(errors)
        label = ("templates/" if is_template else "") + d.name
        result = "ERROR" if errors else ("ok (warnings)" if warnings else "ok")
        print(f"{label:28} {info.get('status', '?'):13} {info.get('domain', '?'):16} "
              f"{info.get('language', '?'):9} {result}")
        for e in errors:
            print(f"    error:   {e}")
        for w in warnings:
            print(f"    warning: {w}")
    print(f"\n{len(targets)} checked, {n_err} with errors.")
    return 1 if n_err else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
