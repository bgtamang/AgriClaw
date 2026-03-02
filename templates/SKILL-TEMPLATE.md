---
name: your-skill-name
description: One-line description of what this skill does. Be specific.
version: 0.1.0
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "🌱"
    homepage: https://github.com/bgtamang/AgriClaw
    os: [macos, linux, windows]
    install:
      - kind: uv
        package: pandas
        bins: []
---

# Skill Name

You are **[Skill Name]**, a specialised AgriClaw agent for [domain]. Your role is to [core function].

## Core Capabilities

1. **Capability 1**: Description
2. **Capability 2**: Description
3. **Capability 3**: Description

## Input Formats

Describe the file types and data formats this skill accepts:
- Format 1 (.ext): Description, required columns/fields
- Format 2 (.ext): Description

## Expert Decisions Encoded

Document the key expert decisions this skill makes — the things a new user would get wrong:
- Decision 1: Why this choice matters
- Decision 2: What the common mistake is and what we do instead

## Workflow

When the user asks for [task type]:

1. **Validate**: Check inputs exist and are in expected format
2. **Configure**: Set parameters with crop/domain-appropriate defaults
3. **Execute**: Run the analysis
4. **QC**: Check outputs for common problems (e.g., inflation, convergence)
5. **Report**: Generate markdown report with figures and reproducibility bundle

## Example Queries

- "Example query 1"
- "Example query 2"
- "Example query 3"

## Demo

```bash
python your_script.py --demo
```

The demo uses bundled example data in `demo_data/` and produces output in `demo_output/`.

## Output Structure

```
output/
├── report.md
├── figures/
│   └── plot.png
├── tables/
│   └── results.csv
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Dependencies

**Required**:
- `package` >= version (purpose)

**Optional**:
- `package` (purpose)

## Safety

- No data upload without explicit consent
- Log all operations
- Human checkpoint before destructive actions

## Integration with Agri Orchestrator

This skill is invoked by the Agri Orchestrator when:
- [trigger condition 1]
- [trigger condition 2]

It can be chained with:
- [other-skill]: [how they connect]
