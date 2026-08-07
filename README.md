# Stomatology-PVL

**A 101-gene direction-weighted progression signature for oral leukoplakia (OL) / proliferative verrucous leukoplakia (PVL): reanalysis and cross-cohort validation of public transcriptomic data.**

This repository contains the full analysis pipeline, results, manuscript draft, and data-deposit package for a reanalysis of public Gene Expression Omnibus (GEO) transcriptomic cohorts. No new experimental data were generated; every result is reproducible from the public datasets and the scripts in this repository.

## Datasets

| Accession | Platform | Samples | Role in this study |
|---|---|---|---|
| [GSE26549](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26549) | Affymetrix Human Gene 1.0 ST (GPL6244) | 86 oral leukoplakia lesions; 62 baseline biopsies (21 progressing / 41 non-progressing) | **Discovery** — 101-gene progression DEG signature (`adj.P < 0.05`, `|logFC| > 0.5`) |
| [GSE227919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227919) | RNA-seq (GPL18573), raw TPM scale | 66 samples: Control 18, HkNR 17, Dysplasia 22, Cancer 9 | **External validation** — severity-gradient correlation of the signature score |
| [GSE196296](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196296) | Single-cell RNA-seq (10x Genomics) | 3 tissue blocks of one PVL lesion (KUS / corrugated / bulky), ~28,000 cells | **Single-cell layer** — cell-state enrichment of the signature |
| [GSE85195](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85195) | Agilent Whole Human Genome 4×44K (GPL6480) | 50 samples: 15 OPL, 34 early OSCC, 1 normal (gingivobuccal complex) | **Second independent cohort** — cross-sectional OPL-vs-OSCC discrimination |

## Key results (honest summary)

- **Signature**: 101 direction-weighted genes derived from GSE26549 progression DEGs; 86 of 101 overlap the previously published risk-model genes.
- **GSE227919 (RNA-seq)**: severity-gradient Spearman ρ = 0.473 (*P* = 6 × 10⁻⁵); ρ = 0.356 after excluding the cancer group (*P* = 0.007) — the signal is not driven by the cancer label alone.
- **GSE85195 (Agilent, different population)**: OPL-vs-OSCC discrimination AUC = 0.949 (95% CI 0.891–1.0), 83/101 genes mapped. **This is a cross-sectional state-discrimination result, not a prospective transformation-prediction result.**
- **Single-cell (GSE196296)**: no lesion-level gradient (between-block η² = 0.006, below the small-effect threshold); the signature is enriched in specific cell states (clusters 5/13/22/23/15, Cliff's δ = 0.43–0.89).
- **Drug reversal (LINCS L1000)**: 8 candidate reversal compounds identified by bidirectional enrichment; PI3K/mTOR inhibitors prioritised by an existing clinical precedent.
- **Caveat**: the internal 5-fold CV-AUC of 1.0 is reported candidly as an optimistic estimate (gene-selection loop bias), not as a generalisation estimate.

## Repository structure

```
PVL_0*.R            Analysis pipeline scripts (R 4.3.3)
*.py                GEO download / metadata parsing scripts
manuscript/         manuscript_v1.md (IMRaD draft) + deposit package
results/            Analysis output tables and figures (per-phase CSVs/PNGs)
docs/               Per-phase result documents (Phase2–8, A1, PVL9)
README.md           This file
```

## Reproducibility

- Each `PVL_0*.R` script reads the raw data and writes to `results/`; the pipeline order is documented in `docs/`.
- The complete file inventory, per-sample scores, and reproduction procedure are in [`manuscript/deposit/README.md`](manuscript/deposit/README.md).
- The deposit package (signature table, per-sample scores for both validation cohorts, single-cell cluster statistics, all scripts, session info) is archived at Zenodo: **https://doi.org/10.5281/zenodo.21823584**

## Citation

If you use this repository or its results, please cite the Zenodo record (doi:10.5281/zenodo.21823584) and the underlying public datasets listed above.

## License

Repository content: [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/). The underlying datasets are subject to their original (GEO) terms.

## AI assistance disclosure

This study is a reanalysis of public data performed with AI-assisted computational tools. Large language models contributed to analysis code, manuscript drafting, and this documentation. All statistical results are generated from the public datasets by the scripts in this repository and were verified by the authors, who take full responsibility for the scientific content. No AI model generated any of the underlying data.
