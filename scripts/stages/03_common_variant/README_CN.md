# MAR03 Common-Variant Retest v1.0R2 REAL

This package supersedes MAR03 R1 before any association-score release.

Why R2 exists:
1. The real-interface probe confirmed the historical MAGMA GSA CLI:
   - marginal: `--model direction=greater`
   - conditional: `--model analyse=file,<target> condition-hide=<others> direction=greater`
2. The probe confirmed the historical MDV8 Match-A1 method:
   chromosome-stratified exponential-tilt rerandomization with exact chromosome composition,
   structural covariates `log_gene_length_z`, `GC_z_mdv8`, and
   `log1p_annotation_degree_z`, per-set |SMD| caliper 0.25, B=10000,
   seed=20260814, and the same matched sets scored in PGC and SPARK.
3. The earlier provisional R1 matched contract (K=50, length+NSNPS) did NOT match
   the frozen historical MDV8 method and therefore must not be released.
4. R2 preserves the already-frozen scientific hypotheses:
   7 higher-order Expanded units, 13 standalones, 7 Compact sensitivities.

Important:
- Do NOT release the old R1 preanalysis lock.
- R2 writes to a clean sibling output:
  `/home/h3021/chapter3/12_molecular_autism_revision/03_common_retest/MAR03_common_variant_retest_v1.0R2`
- No SNP harmonization or MAGMA gene mapping is rerun.
- Existing PGC/SPARK gene-level authorities are reused.
- No `conda run` is used. The runner calls the already validated Python executable directly.

## Execution

```bash
cd /home/h3021/chapter3/12_molecular_autism_revision/MAR03_COMMON_VARIANT_RETEST_v1.0R2_REAL

bash run_mar03_fast.sh all
```

Expected:
`MAR03_PREFLIGHT=PASS`
`MAR03_FREEZE=PASS`
`A_B_C=7/13/7`
`HOLD_FOR_MANUAL_GATE`

Then inspect:
- `MAR03_preflight_report.tsv`
- `MAR03_preanalysis_coverage.tsv`
- `MAR03_preanalysis_lock.json`
- `MAR03_HISTORICAL_INTERFACE_SNAPSHOT.tsv`

If correct:
```bash
bash run_mar03_fast.sh release RELEASE
bash run_mar03_fast.sh analysis
bash run_mar03_fast.sh matched
bash run_mar03_fast.sh sensitivity
bash run_mar03_fast.sh postflight
bash run_mar03_fast.sh pack-review
```

Do not re-run `all` after release because freeze timestamps would change the lock hash.

## Frozen scientific families

- Family A: 7 higher-order Expanded units, primary competitive MAGMA, BH m=7.
- Family B: 13 standalones, secondary completeness, BH m=13.
- Family C: 7 Compact units, sensitivity.

Conditional A-unit tests are supportive only and remain
`SUPPORTIVE_PENDING_MAR05`.

## Matched continuous analysis

The matched analysis is score-blind during set construction.
All frozen MAR03 A+B hypothesis genes are excluded from the control pool.
The same 10,000 unique matched sets per target are scored in PGC and SPARK.

Mean Z:
- primary complementary statistic.
- empirical P uses `(1 + #null >= observed)/(B+1)`.
- BH is applied separately within dataset × family.

Median Z:
- supportive.
- same empirical and within-family BH treatment.

## Sensitivities

For the 7 higher-order units:
- Compact w10
- Expanded remove-raw26 w10
- Expanded w0
- Expanded w50
- Expanded MHC-excluded

Expected sensitivity rows: 2 datasets × 7 units × 5 scenarios = 70.
