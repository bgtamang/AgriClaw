# Test fixtures

- `mini.gff3`: hand-made GFF3 (invented genome "AgriClaw-mini-1.0") with genes placed on the boundaries of the interval Gm01:1000-2000: ending 1 bp before it, touching its start, inside it, starting at its end, starting 1 bp after it and spanning it. It also has a pseudogene, a gene on another chromosome, URL-encoded attributes (`%2C`, `%3B`), `gene:` ID prefixes as in Ensembl GFF3 files, and a feature line after `##FASTA` that must be ignored.
- `ensembl_info_assembly_glycine_max.json` and `ensembl_overlap_glycine_max_1_1000000-1060000.json`: responses saved from the Ensembl REST API (https://rest.ensembl.org, Ensembl release 116) on 2026-10-15:
  - `GET /info/assembly/glycine_max`, with `top_level_region` trimmed to the 20 chromosomes to keep the file small;
  - `GET /overlap/region/glycine_max/1:1000000-1060000?feature=gene`.

  The soybean gene models are the JGI Wm82.a2.v1 annotation on assembly Glycine_max_v2.1, as served by Ensembl Plants. Ensembl states that its data are available without restriction (https://www.ensembl.org/info/about/legal/disclaimer.html); the responses are used here only to test JSON parsing offline.
