# Data and Code Deposit — Direction-weighted 101-gene progression signature for oral leukoplakia / PVL

**Version 2** · Zenodo record 21823584 · concept DOI 10.5281/zenodo.21823583

This package contains the derived data tables, per-sample scores and analysis code required to
reproduce the results of the manuscript. The study is a re-analysis of four public GEO datasets
(GSE26549, GSE227919, GSE196296, GSE85195); no new measurements were generated.

---

## What changed in version 2

Version 1 was released on 2026-08-06. Version 2 fixes defects that made version 1 **not
reproducible as documented**, and adds the source tables behind every reported number.

| Change | Detail |
|---|---|
| **Scripts updated** | All nine `PVL_0*.R` scripts were the pre-relocation versions. They hard-coded `/mnt/g/...` paths and did not read the configuration file at all, so the `DATASET_PATH` variable that the version-1 README told readers to set was never actually consulted. They are replaced by the current versions. |
| **`_env_paths.R` added** | The single source of truth for paths, sourced by every analysis script. It was missing from version 1. |
| **`.env.template` fixed** | Version 1 contained only comment lines and defined no variables. It now carries real `DATASET_PATH` / `DATASET_PATH_TEMP` entries. |
| **`signature_101genes.csv` corrected** | The `signif` column had been computed with a stricter `|logFC| > 1.0` boundary, so only 13 of 101 rows read `yes`. It is now computed with the documented threshold (`adj.P < 0.05` and `|logFC| > 0.5`) and all 101 rows read `yes`. |
| **`PVL6_summary.csv` corrected** | `KW_p` had been formatted into numeric underflow and stored as `0`. It now reads `1.127e-11`, with `KW_chisq` (50.418) and `KW_df` (2) added so the value can be checked. |
| **Source tables added** | Version 1 omitted the tables behind Table 1, the AUCs, the probe-collapse sensitivity analysis and the cluster statistics. All are now included (see inventory). |
| **Single-cell audit tables added** | Cell-type-by-block composition and marker-defined cell-state detection, which let the manuscript's cell-state statements be checked directly. |
| **Audit scripts added** | `PVL_11` (cell-type-by-block composition) and `PVL_12` (marker-defined state check). |
| **`Manifest.csv` added** | SHA-256 checksums of every file in this version. |

---

## File inventory

### Derived data tables (root)

| File | Description |
|---|---|
| `signature_101genes.csv` | The 101-gene progression signature: `gene`, discovery `logFC`, `adj.P.Val`, `AveExpr`, `signif`. Derived from GSE26549 baseline biopsies (21 progressing vs 41 non-progressing) by limma at `adj.P < 0.05` and `|logFC| > 0.5`; all 101 rows satisfy that threshold. |
| `GSE26549_progression_DEG.csv` | Full discovery differential-expression table (all genes), the source of the signature. |
| `per_sample_scores_GSE227919.csv` | Per-sample direction-weighted score for the 66-sample RNA-seq cohort: `sample`, `group` (Control / HkNR / Dysplasia / Cancer), ordinal `severity` (0–3), `direction_weighted_score`. |
| `per_sample_scores_GSE85195.csv` | Per-sample direction-weighted score for the 50-sample Agilent cohort: `sample`, `gsm`, `group` (OPL / OSCC / Normal), `title`, `direction_weighted_score`. |
| `PVL9_GSE85195_scores.csv` | Scores underlying Figure 4 (`sample`, `gsm`, `group`, `title`, `score`). Same values as `per_sample_scores_GSE85195.csv` with the column named `score`, kept because the figure script reads this file. |
| `PVL9_summary.csv` | Second-cohort summary metrics: sample counts, mapped gene count (83/101), AUC (0.949), DeLong CI, severity ρ. |
| `PVL9_mapped_genes.csv` | The 83 signature genes that map to the GSE85195 features, with the contributing Agilent probes. Rebuilds the mapping that version 1 reported only as a count. |
| `PVL4_summary.csv` | Risk-model metrics: internal 5-fold CV AUC, direct LASSO transfer AUC, progression-subset AUC, pseudobulk direction concordance. |
| `PVL4_lasso_coefficients.csv` | The 20 LASSO-selected genes and their coefficients. |
| `PVL5_signature_signed.csv` | The signature split by direction as submitted to LINCS: `up_in_progression` (26 genes) and `down_in_progression` (60 genes). **86 of the 101 genes are mappable to L1000; the remaining 15 are multi-symbol / lncRNA / LOC entries with no single L1000 gene symbol.** |
| `PVL5_lincs_reversal_candidates.csv` | The eight two-direction reversal candidates with per-direction scores, overlapping genes and the consensus-layer flag. |
| `PVL6_summary.csv` | Single-cell summary: per-block mean score, Kruskal–Wallis statistic (`KW_chisq`, `KW_df`, `KW_p`), between-block η², top cluster. |
| `PVL6_lasso_marker_overlap.csv` | The four LASSO genes that overlap epithelial/keratinocyte cluster markers. |
| `PVL8_probe_sensitivity.csv` | Probe-collapse sensitivity: signature size under max vs mean collapsing, Jaccard, logFC Spearman ρ, external ρ. |
| `PVL10_cluster_sigscore_stats.csv` | Cluster-level localisation statistics: n cells, mean score, Wilcoxon p, Cliff's δ, η², FDR for all 24 clusters. |
| `PVL10_summary.csv` | Cluster and cell totals. |
| `SC_celltype_by_lesion_counts.csv` | Cell-type × lesion-block cell counts (24 cell types). |
| `SC_celltype_by_lesion_props.csv` | The same table column-normalised to percentages within each block. |
| `SC_celltype_by_lesion_stats.csv` | Per-cell-type test across blocks: proportions, `prop.test` p, Cramér's V, FDR. |
| `SC_marker_detection_by_lesion.csv` | Per-block detection rate and mean expression for 32 marker genes used to define immune cell states. |
| `SC_marker_state_by_lesion_stats.csv` | Marker-defined cell states (Treg, plasma cell, mast cell) per block with p, Cramér's V and FDR. |
| `session_info.txt` | `sessionInfo()` of the R environment that produced the results (R 4.3.3, Seurat v5, package versions). |
| `Manifest.csv` | SHA-256 of every file in this version. |

### Code (shipped alongside the data tables, flat layout)

| File | Description |
|---|---|
| `_env_paths.R` | Reads `DATASET_PATH` and `DATASET_PATH_TEMP` from the project `.env` and exports `RAW`, `MET`, `RES`. Sourced by every analysis script. |
| `env.template` | Template for the `.env` file. Copy it to `.env` in the same directory and edit the two paths. |
| `PVL_02_bulk.R` … `PVL_10_sc_cluster_stats.R` | The analysis pipeline, in execution order (see below). |
| `PVL_11_sc_lesion_celltype_audit.R` | Cell-type-by-block composition and testing. |
| `PVL_12_sc_marker_claim_check.R` | Marker-defined cell-state detection across blocks. |
| `R_install.R`, `R_install2.R` … `R_install5.R` | Package installation scripts. |
| `download_*.py`, `fetch_metadata.py`, `parse_*.py` | Utilities that fetched the GEO series matrices and parsed the clinical tables. |

---

## How to reproduce

### 1. Obtain the public input data

Download from GEO into `$DATASET_PATH`:

| Dataset | Platform | Put in |
|---|---|---|
| GSE26549 | Affymetrix Human Gene 1.0 ST (GPL6244) | `raw/GSE26549_series_matrix.txt.gz` plus annotation `raw/GPL6244.annot.gz` |
| GSE227919 | RNA-seq (GPL18573) | `raw/GSE227919_*.txt.gz`, `raw/GSE227919_series_matrix.txt.gz`, `raw/GSE227919_2023_21_03_pml_annotations.csv.gz` |
| GSE196296 | 10x scRNA-seq | `raw/GSM5870060_LesionA_*.gz`, `raw/GSM5870062_LesionB_*.gz`, `raw/GSM5870063_LesionC_*.gz` |
| GSE85195 | Agilent 4×44K (GPL6480) | `raw/GSE85195/GSE85195_series_matrix.txt` |

The scripts `download_full.py`, `download_sc.py`, `download_suppl.py` and `fetch_metadata.py`
automate this step. Parsed clinical tables (`GSE26549_clinical.csv`, `GSE227919_clinical.csv`,
`GSE196296_sampletable.csv`) are written to `metadata/`.

### 2. Configure the paths

```bash
cp env.template .env      # then edit DATASET_PATH and DATASET_PATH_TEMP
```

Both values are Windows-style paths (e.g. `E:\Medical Data\Stomatology\PVL`).
`_env_paths.R` translates them to WSL mounts (`/mnt/e/...`) when run inside WSL.
Run the analysis scripts from the **project root**, because they locate `_env_paths.R` and
`.env` via `getwd()`. If no `.env` is present the scripts stop with an explanatory message.

### 3. Install the R environment

```bash
Rscript R_install.R
Rscript R_install2.R
Rscript R_install3.R
Rscript R_install4.R
Rscript R_install5.R
```

### 4. Run the pipeline, in this order

```bash
Rscript PVL_02_bulk.R               # discovery DEGs and the 101-gene signature
Rscript PVL_03_scrna.R              # single-cell object construction
Rscript PVL_04_riskmodel.R          # LASSO risk model, cross-validation, external transfer
Rscript PVL_05_lincs.R              # directional LINCS queries via enrichR
Rscript PVL_06_closeloop.R          # single-cell localisation of the signature
Rscript PVL_07_signature_decon.R    # signature-level enrichment
Rscript PVL_08_probe_sens.R         # probe-collapse sensitivity
Rscript PVL_09_second_cohort.R      # GSE85195 second cohort
Rscript PVL_10_sc_cluster_stats.R   # cluster-level effect sizes
```

`PVL_11` and `PVL_12` are audit scripts and can be run after `PVL_10`. They require the
single-cell object from `PVL_03` (`$DATASET_PATH_TEMP/PVL_sc_combined.rds`).

Every script writes to `$DATASET_PATH_TEMP`. It does not write to the package directory.

---

## Notes on specific numbers

- **Signature size.** 101 genes are derived from GSE26549. 86 of them are present in GSE227919.
  Of the 101, 86 are mappable to L1000 and those are the genes submitted to the LINCS queries
  (26 progression-up, 60 progression-down); the other 15 have no single L1000 gene symbol.
- **LINCS candidate count.** The two directional queries at `adj.P < 0.05` return 720 terms
  (down-gene set against `Chem_Pert_up`) and 8 terms (up-gene set against `Chem_Pert_down`).
  The two sets share exactly one term, so the union is 727 term entries, corresponding to
  244 unique compounds. Requiring significance in both directions leaves eight compounds.
- **Consensus-layer flag.** `consensus_confirmed` in `PVL5_lincs_reversal_candidates.csv` is
  `FALSE` for all eight. This is expected and is not evidence against the candidates: the
  `LINCS_L1000_Chem_Pert_Consensus_Sigs` library currently exposes 1,413 terms, all in the `Up`
  direction and none in `Down`, while the protective direction requires `Down`. Seven of the
  eight compounds are not represented in that library at all.
- **Cross-validation AUC.** The internal five-fold cross-validated AUC of 1.0 is optimistic,
  because gene selection and evaluation are not nested. It is reported as an upper bound, not
  as a performance estimate. The independent evidence is the GSE227919 severity gradient and the
  GSE85195 cross-sectional discrimination.
- **GSE85195 result.** The OPL-versus-OSCC discrimination (AUC 0.949) is a cross-sectional
  state-discrimination result on a different population (gingivobuccal complex). It is not a
  prospective transformation prediction.
- **Single-cell gradient.** The between-block η² of 0.006 is below the small-effect threshold
  used in the manuscript, so the signature does not track the KUS-to-bulky axis. The
  Kruskal–Wallis p value is reported alongside the effect size precisely to make that point.
- **Enrichr library drift.** The directional LINCS libraries are live. Term counts can change by
  a few entries between Enrichr releases; the numbers above were obtained in October 2026.

---

## Licence

- Analysis code (all R and Python files): MIT.
- Data tables (`*.csv`, `session_info.txt`): CC-BY-4.0.

The four input datasets remain under their original GEO/NCBI terms.

## Citation

If you use this package, please cite the Zenodo record (DOI 10.5281/zenodo.21823584, or the
concept DOI 10.5281/zenodo.21823583 which always resolves to the latest version) and the
underlying public datasets listed above.
