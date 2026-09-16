# IEI-ASD NeuroGenomics

Publication reproducibility repository for a Molecular Autism study examining how far pathway hypotheses derived from an inborn-errors-of-immunity (IEI)–autism gene-catalog intersection remain supportable across pathway prioritization, statistical calibration, and common-variant analyses.

## Scientific scope

A curated 26-gene IEI–autism catalog intersection was used as a hypothesis-generating anchor rather than as evidence of patient-level rare-variant convergence.

The rare-side analyses yielded:

- 18 pathways passing the dual-background comparison;
- 39 pathways prioritized under the prespecified BH-based analysis;
- 30 nonredundant representatives;
- 7 higher-order common-variant hypotheses;
- 13 secondary standalone hypotheses.

The seven higher-order hypotheses showed no reproducible preferential common-variant enrichment in either PGC2019 or SPARK.

A dependency-robust multiplicity sensitivity retained none of the 39 pathway priorities. Full-chain global-null replay further showed a high probability of at least one false pathway-priority event under the complete null.

Accordingly, the rare-side outputs are treated as a descriptive, redundancy-compressible hypothesis space rather than as validated mechanisms.

## Interpretation boundary

This repository does **not** support:

- a single unified IEI–autism neuroimmune mechanism;
- patient-level co-occurrence of the catalog variants;
- absence or irrelevance of common-variant effects;
- clinical risk stratification;
- treatment selection or biomarker deployment.

The study instead defines where catalog-derived pathway hypotheses remain supportable and where stronger interpretation stops.

## Repository organization

- `protocol/` — analysis locks, contracts, and provenance records.
- `scripts/` — final analysis code used for the reported results.
- `source_data/` — publication-level derived data used for figures, tables, and reported numerical summaries.
- `figures/` — final vector manuscript figures.
- `environment/` — software/environment documentation.
- `tests/` — automated checks for key manuscript claims.
- `checksums/` — integrity manifests.
- `docs/` — study overview, provenance, portability, and claim boundaries.

## Archival release

The Molecular Autism submission snapshot is archived on Zenodo:

**Version 1.0.0:** DOI: 10.5281/zenodo.22790971


## External data

Primary third-party genomic resources are not redistributed here. These include PGC and SPARK association data and the 1000 Genomes reference resource. Users should obtain these data from the original providers under their applicable access and redistribution conditions.

Derived publication-level summary tables required to verify the manuscript results are included in `source_data/`.

## Reproducibility

The repository preserves the exact analysis-code snapshots and analysis-lock records used for the manuscript. Some scripts retain historical HPC paths as fallback defaults. Portable execution is supported through environment variables such as `CH3_ROOT` and `MA_ROOT`; see `docs/PORTABILITY.md`.

The continuous-integration workflow validates:

- Python syntax;
- native R parse gates;
- shell syntax;
- key manuscript claim counts;
- repository checksums.

## Version history

The earlier exploratory CoreSeed/three-axis repository state is preserved in:

- branch: `legacy-coreseed-axis`
- tag: `legacy-pre-MDV-MAR-20260915`

Those materials are historical and are not the analytic authority for the current manuscript.

## License

The MIT license applies to repository code. External datasets and derived materials remain subject to the terms of their original data providers.
