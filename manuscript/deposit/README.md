# Data and Code Deposit — "A validated transcriptomic signature of malignant progression in oral leukoplakia and PVL"

This package contains all data objects and analysis code required to reproduce the
results of the manuscript. The study is a re-analysis of four public GEO datasets
(GSE26549, GSE227919, GSE196296, GSE85195); no new measurements were generated.

## File inventory

| File | Description |
|------|-------------|
| `signature_101genes.csv` | The 101-gene progression signature (gene, discovery logFC, adj.P.Val, AveExpr, signif). Derived from GSE26549 (limma; adj.P<0.05 & \|logFC\|>0.5). |
| `GSE26549_progression_DEG.csv` | Full discovery differential-expression table (all genes), source of the signature. |
| `per_sample_scores_GSE227919.csv` | Per-sample direction-weighted z-score for the 66-sample RNA-seq validation cohort (GSE227919): sample, group (Control/HkNR/Dysplasia/Cancer), ordinal severity (0–3), score. |
| `per_sample_scores_GSE85195.csv` | Per-sample direction-weighted z-score for the 50-sample Agilent cross-sectional cohort (GSE85195): sample, gsm, group (OPL/OSCC), title, score. |
| `session_info.txt` | `sessionInfo()` of the R environment that produced the results (R 4.3.3, Seurat v5, package versions). |
| `scripts/` | All analysis scripts (PVL_02 … PVL_10) and environment install scripts (R_install*.R), plus the download/parse Python utilities and a sanitised `.env.template`. |
| `upload_to_zenodo.py` | Helper that publishes this package to a Zenodo repository via the REST API (requires a Zenodo access token). |

## How to reproduce

1. Set `DATASET_PATH` in `scripts/.env` (copy from `.env.template`) to the directory
   that contains the downloaded `raw/` and `metadata/` subdirectories of the four GEO
   datasets. On WSL, this is translated to `/mnt/<drive>/...`.
2. Install R packages: run `scripts/R_install.R` (and `R_install2.R` … `R_install5.R`).
3. Run the pipeline in order: `PVL_02_bulk.R` → `PVL_03_scrna.R` → `PVL_04_riskmodel.R`
   → `PVL_05_lincs.R` → `PVL_06_closeloop.R` → `PVL_07_signature_decon.R`
   → `PVL_08_probe_sens.R` → `PVL_09_second_cohort.R` → `PVL_10_sc_cluster_stats.R`.
   Each script reads from `DATASET_PATH` and writes CSV/PNG summaries to its results dir.

## Licence

- Analysis code (`scripts/`): MIT.
- Data tables (`*.csv`): CC-BY-4.0.

## Publishing to Zenodo

```bash
export ZENODO_TOKEN=your_personal_access_token
python3 upload_to_zenodo.py --sandbox   # dry-run against sandbox first
python3 upload_to_zenodo.py             # publish to production; prints the DOI
```

The deposit DOI is https://doi.org/10.5281/zenodo.21823584 (inserted into the
manuscript Data Availability Statement). Confirm the record is in "Published" state
on Zenodo — a reserved DOI returns HTTP 404 from doi.org until the record is published.
