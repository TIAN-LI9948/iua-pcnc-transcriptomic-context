# Build provenance

- Build date: 2026-10-31 (Asia/Shanghai)
- Primary statistical source: audited project script `22_phase_C_D_pc_nc_v2.R`
- Cross-context source: audited project script `24_phase_E_cross_context_v2.R`
- Untimed bulk boundary source: audited project script `31_phase_I_GSE224093_bulk_contract_and_endpoints.py`
- Paired-extension sources: `12_reaggregate_gse311899_paired_epithelial.py` and `13_test_gse311899_paired_longitudinal.R`
- Portable reproduction entry point: `scripts/01_reproduce_primary_endpoints_figure1_table2.R`
- Numerical validation entry point: `scripts/02_validate_reproduction.py`
- R runtime used for the validated rerun: R 4.6.0 (see `reproduced_outputs/R_sessionInfo.txt`)
- Python runtime used for validation: Python 3.11.13 (see `metadata/reproduction_validation.json`)
- Input/output fingerprints: `metadata/file_manifest_sha256.tsv`
- Numerical QA: `metadata/reproduction_validation.json` (`PASS`, 20/20 endpoint-field comparisons within `1e-12`)
- Runtime/package versions: `reproduced_outputs/R_sessionInfo.txt` and the Python runtime field in `metadata/reproduction_validation.json`
- Upstream scripts that require the complete project workspace read its location from the `IUA_CONTEXT_ROOT` environment variable; the portable primary reproduction does not require this variable.

The portable script was checked against the frozen primary result tables. All four endpoint differences, bootstrap confidence limits, exact permutation *P* values, and within-family Benjamini–Hochberg *q* values matched exactly at the recorded tolerance.
