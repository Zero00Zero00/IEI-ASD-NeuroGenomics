# MAR07_CHAIN_CALIBRATION_v1.0R4.1_NA_SAFE_CERTIFICATION

R4.1 fixes only the replay-comparator implementation bug observed at `certify`:
the R MDV2 replay itself completed successfully, but the Python comparator attempted `float("NA")`.

Scientific design is unchanged. For reproducibility R4.1 uses a fresh output namespace:
`07_chain_calibration/MAR07_chain_calibration_v1.0R4p1`
and therefore requires a fresh freeze/release.

Run:
1. `bash run_mar07_r4p1.sh all`
2. inspect `MAR07_preanalysis_lock.json`
3. `bash run_mar07_r4p1.sh release RELEASE`
4. `bash run_mar07_r4p1.sh observed`
5. `bash run_mar07_r4p1.sh certify`
6. only after `MAR07_FULLCHAIN_EXECUTOR_CERTIFIED=PASS`:
   `bash run_mar07_r4p1.sh pilot`
7. `bash run_mar07_r4p1.sh calibrate`
8. `bash run_mar07_r4p1.sh postflight`
9. `bash run_mar07_r4p1.sh pack-review`

New diagnostic:
`certification/MAR07_CERT_NUMERIC_NA_AUDIT.tsv`

No solver/model/FDR/threshold/pathway-family change was made.
