# Data provenance

## Principle

This repository separates third-party primary genomic resources from publication-level derived data.

The repository redistributes only code, analysis metadata, and derived source tables that are required to verify the reported manuscript results and that are appropriate for public redistribution.

## Catalog and annotation resources

The analysis was constructed from curated IEI and high-confidence autism gene catalogs mapped into a common protein-coding gene universe.

Pathway analyses used Gene Ontology Biological Process and Reactome membership after the project-specific filtering and harmonization steps documented in the analysis code and lock files.

Exact input filenames, hashes, mapping decisions, pathway filters, and gene-universe definitions are preserved in the protocol and checksum records.

## Derived inputs distributed in this repository

`source_data/derived_inputs/` contains publication-relevant derived inputs including:

- the protein-coding analysis universe;
- derived gene-group assignments;
- pathway manifest information;
- the dual-background pathway candidate table.

These files are included to permit audit of manuscript-level results.

## Association resources not redistributed

The following primary or high-dimensional resources are not redistributed:

- PGC2019 SNP-level autism GWAS summary statistics;
- SPARK SNP-level autism GWAS summary statistics;
- full harmonized SNP-level intermediate data;
- 1000 Genomes genotype/LD reference data;
- full per-gene association outputs where redistribution conditions have not been independently confirmed.

Users must obtain these resources from the original providers under their applicable terms.

## Common-variant publication data

The repository distributes only publication-level aggregate quantities required to verify the manuscript, including:

- gene-set beta estimates;
- standard errors;
- one-sided P values;
- BH-adjusted q values;
- coverage summaries;
- matched-continuous empirical percentiles;
- precision and MDE summaries;
- collinearity diagnostics.

## Global-null calibration

The full high-dimensional replicate working directories are not distributed.

The repository includes the replicate-level stage counts and calibration summaries required to reproduce the reported quantities, including:

- the number of priority events per replicate;
- the probability of at least one false priority event;
- the mean null priority count;
- Monte Carlo uncertainty summaries.

## Driver-influence analyses

Publication-level leave-one-driver and leave-two-driver summaries are distributed. These are diagnostic sensitivity analyses and do not redefine the upstream pathway-prioritization classifications.

## Historical absolute paths

Some archived analysis scripts retain original `/home/h3021/...` paths as fallback defaults. They are preserved to maintain identity with the executed code. Portable execution should use environment variables and configuration files as documented in `PORTABILITY.md`.
