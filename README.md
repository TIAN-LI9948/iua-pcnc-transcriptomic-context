# Reproducibility materials for the IUA/PC–NC transcriptomic-context study

This repository supports the secondary public-data analysis reported in:

> *Epithelial ESR1/PGR Expression and a Window-of-Implantation Program in Programmed-Cycle and Intrauterine Adhesion Contexts: A Secondary Analysis of Public Transcriptomic Data*

The primary analysis compares donor-level nonciliated epithelial pseudobulks from 7 programmed-cycle (PC) and 9 natural-cycle (NC) participants with unexplained infertility. The Asherman syndrome/hormone replacement therapy (AS/HRT) dataset is a composite contextual comparison: disease status and regimen are aligned, so these materials must not be used to claim a disease-only or HRT-adjusted effect. The three eligible pre/post GSE311899 pairs are descriptive and overlap the source population of the AS study; they are not an independent validation or a treatment-effect analysis.

## Quick reproduction of Figure 1 and Table 2

Required software:

- R 4.6.0 or a compatible recent R release
- R packages: `edgeR`, `data.table`, `ggplot2`, and `cowplot`
- Python 3.9 or later for the dependency-free validation script

From the repository root, run:

```powershell
Rscript scripts/01_reproduce_primary_endpoints_figure1_table2.R reproduced_outputs
python scripts/02_validate_reproduction.py
```

On Windows, an explicit path to `Rscript.exe` may be used if R is not on `PATH`.

The R script regenerates:

- `reproduced_outputs/primary_endpoint_scores.tsv`
- `reproduced_outputs/primary_endpoint_tests.tsv`
- `reproduced_outputs/Table_2_endpoint_estimates.tsv`
- `reproduced_outputs/Figure_1.pdf`
- `reproduced_outputs/Figure_1.png`
- `reproduced_outputs/R_sessionInfo.txt`

The validation script compares every primary numerical estimate, bootstrap interval, exact permutation *P* value, and within-family Benjamini–Hochberg *q* value with the frozen expected outputs at an absolute tolerance of `1e-12`. It writes `metadata/reproduction_validation.json` and a SHA-256 file manifest.

## Primary scoring contract

- Biological inference unit: donor, not nucleus.
- Compartment: nonciliated epithelium (publisher glandular and luminal epithelial labels; ciliated epithelial and cervical epithelial labels excluded).
- Normalization: edgeR TMMwsp library normalization followed by logCPM with prior count 2.
- PGR and ESR1: donor-pseudobulk transcript abundance in logCPM.
- Frozen WOI program: unweighted arithmetic mean of 22 mapped genes from a 25-gene program defined outside the PC/NC target dataset. No target-dataset gene selection, gene-level standardization, or imputation is performed.
- Stable epithelial core: unweighted arithmetic mean of 9 mapped genes from an 11-gene reference stability core. No target-dataset gene selection, gene-level standardization, or imputation is performed.
- Exact test: all 11,440 allocations of 7 PC labels among 16 donors.
- Confidence interval: 20,000 within-group nonparametric donor-bootstrap replicates.
- Multiplicity: PGR plus frozen WOI program form the primary family; ESR1 plus stable epithelial core form the secondary family.

Complete program membership, mapping, direction, and effective weights are provided in `data/derived/score_gene_mapping_and_direction.tsv`.

## Directory guide

- `data/primary/`: donor-pseudobulk integer counts, donor metadata, and the two externally defined program contracts needed for portable primary reproduction.
- `data/derived/`: contextual endpoint estimates, donor metadata/missingness, design-sensitivity results, and the gene-level scoring contract.
- `data/figure_source/`: source tables used by the final main and supplementary figures.
- `expected_outputs/`: frozen primary results used only for numerical validation.
- `scripts/`: portable primary reproduction and validation entry points.
- `scripts/validated_pipeline/`: preserved analysis and plotting scripts from the audited project pipeline. These document the upstream cohort contracts, cross-context comparison, bulk boundary analysis, paired extension, and final figure rendering. Upstream scripts that require the complete raw-data workspace read its location from the `IUA_CONTEXT_ROOT` environment variable. Use the portable script above for direct reproduction from the included inputs.
- `metadata/`: checksums and machine-readable reproduction QA.

## Public data sources

The original inputs are publicly accessible from CELLxGENE Discover collection `e80f5450-0a0f-4312-a728-431d50aeb8fe` and NCBI Gene Expression Omnibus accessions `GSE111976`, `GSE200805`, `GSE215968`, `GSE224093`, `GSE250130`, and `GSE311899`.

The repository contains derived donor-level or figure-source tables needed to audit the reported results; it does not redistribute the large original single-cell objects. Original accessions, cohort roles, inference units, and interpretation limits are reported in the manuscript and supplementary tables.

## Evidence limits

These scripts reproduce molecular associations at the sampled time points. They do not measure receptor protein abundance or activity, endometrial receptivity, implantation, pregnancy, live birth, or treatment efficacy. The opposite stable-core direction in the untimed GSE224093 bulk dataset and the nonsignificant branches are retained rather than filtered out.

## License and reuse

Project-authored code is released under the MIT License. Derived tables are provided for transparent verification of the reported analyses. The original public datasets remain subject to the terms of their source repositories and studies; users should cite the corresponding source publications and accessions.
