# v2.1 additions — files missing from the published Zenodo v2

**Status: prepared, NOT published.** Published record remains v2 (10.5281/zenodo.23453360).

## Why this directory exists

Zenodo v2 contains 47 payload files: the derived data tables, the full pipeline, the
environment scripts and a SHA-256 manifest. It does **not** contain:

1. the six figure/table source tables under `results/manuscript_figs/`, and
2. the twelve raw Enrichr result files that the bulk-level enrichment numbers derive from.

Consequence: a reviewer cannot verify Figure 3 or the bulk-level columns of Table 2 from
the archive alone — they would have to re-run `PVL_07_signature_decon.R`. That is a
reproducibility gap, not a correctness problem, because the pipeline reproduces the
numbers exactly (verified below).

## What is here

| File | Lines | Backs |
|---|---|---|
| `Fig1A_data.csv` | 66 | Figure 1A (per-sample score, `sample/Class/severity/score`) |
| `Fig2_reversal_matrix.csv` | 8 | Figure 2 (8 compounds × 8 hub genes) |
| `Fig3_data.csv` | 10 | Figure 3 (top Hallmark terms; see the caveat below) |
| `Table1_reversal_candidates.csv` | 8 | Table 1 (condensed form of `PVL5_lincs_reversal_candidates.csv`) |
| `Table2_enrichment_summary.csv` | 2 | Table 2 (signature level vs bulk level) |
| `Table2_bulk_enrichment_detail.csv` | 6 | Table 2 bulk column, per library and direction |
| `enrich_*_GSE227919_CancerVsControl_{up,down}.csv` | 50–5,116 | Figure 3, Table 2 bulk level |
| `enrich_*_GSE26549_progression_{up,down}.csv` | 11–1,569 | Table 2 signature level |

Twelve `enrich_*` files: three libraries (`MSigDB_Hallmark_2020`,
`GO_Biological_Process_2021`, `KEGG_2021_Human`) × two contrasts × two directions.

**Total 18 files, 1.4 MB.**

## Verification performed on this set

| Check | Result |
|---|---|
| Hallmark up, adj.P < 0.05 → count | **24**, matching the manuscript's Table 2 |
| `Table2_bulk_enrichment_detail.csv`, all six rows recomputed from the raw files | 24 / 3 / 405 / 3 / 59 / 0 — every row matches |
| Rebuilding `Fig3_data.csv` top terms from the Hallmark raw file | first 9 terms match by combined score; **10th differs** — see caveat |
| Hard-coded host paths in any script | none |

## Caveats to carry into the manuscript (do not silently archive without noting)

**1. Figure 3 legend vs figure content.** The legend says "top ten significant terms …
ranked by combined score", but the 10th item plotted is `PI3K/AKT/mTOR Signaling`, whose
true rank is **24 of 24** (adjusted P = 0.041). The genuine 10th by combined score is
`Myc Targets V2` (55.96, adjusted P = 4.49e-6). PI3K/AKT/mTOR was added because the text
emphasises it; the legend should say so, or the panel should plot the true top ten. This
is documented in `docs/写作团队交接文案.md` §一 E.

**2. Two kinds of zero in Table 2 — do not conflate.**
- Signature level, TF libraries: **queried and nothing reached significance**
  (`PVL7_ChEA_TF_top.csv` and `PVL7_ENCODE_TF_top.csv` are header-only). Real negative.
- Bulk level, TF libraries: **never queried.** The bulk enrichment step ran only
  Hallmark, GO and KEGG. This is "not analysed", not "not detected".

## How to publish v2.1, if that is decided

The addition would go on top of the existing v2 file set. The DOI workflow is:

1. `POST /api/deposit/depositions/23453360/actions/newversion` → returns the new draft id
2. keep the 47 v2 files (they are inherited automatically), add these 18 via
   `PUT {bucket}/{filename}`
3. regenerate `Manifest.csv` over all 65 payload files
4. amend the README's inventory table and the "What changed" table to describe v2.1
5. `POST /api/deposit/depositions/{draft}/actions/publish`

Environment notes for the Zenodo API: `zenodo.org` is DNS-poisoned on this machine
(resolves to 0.0.0.0) — use `--resolve zenodo.org:443:137.138.121.78 --noproxy '*'`;
file deletion requires the per-file `links.self` UUID, not the filename.

If v2.1 is **not** published, the manuscript's Data Availability statement stays correct
as written (it cites the concept DOI), but the reproducibility gap above remains and
should be acknowledged if a reviewer raises it.

## Decision point

Publishing costs one new version number and a small README edit; not publishing costs a
reviewer round-trip re-running `PVL_07`. Both are defensible. This is the author's call.
