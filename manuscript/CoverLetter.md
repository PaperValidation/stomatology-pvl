2026-10-24

The Editor-in-Chief
Oral Oncology
Elsevier

Dear Editor-in-Chief,

We submit our original research article, *A validated transcriptomic signature of malignant progression in oral leukoplakia and proliferative verrucous leukoplakia: independent replication, single-cell mapping, and in silico drug repositioning*, as an Original Research Article.

Malignant transformation in oral leukoplakia (OL) and proliferative verrucous leukoplakia (PVL) is clinically unpredictable: dysplasia grading is subjective and weakly predictive, and resecting OL has not been shown to reduce the incidence of oral squamous cell carcinoma (OSCC). Saintigny's 29-transcript model and Foy's immunological and classical subtypes, the best-known solutions, both derive from GSE26549, and neither has been independently confirmed on a modern RNA-seq platform. Rather than transferring fitted coefficients across platforms, we score samples with a direction-weighted, within-cohort z-score that borrows only each gene's discovery direction and relative magnitude, then test it in an independent RNA-seq cohort and a second, independent microarray cohort.

Cross-platform validation: a 101-gene signature from GSE26549 (86 OL; 35 progressors, 51 non-progressors) tracked severity in the independent RNA-seq cohort GSE227919 (n = 66; Control 18, HkNR 17, dysplasia 22, cancer 9; ρ = 0.473, P = 5.99 × 10⁻⁵) and survived removal of the nine cancers (n = 57; ρ = 0.356, P = 0.007); in GSE85195 (15 OPL, 34 OSCC; 83/101 genes mapped) it separated pre-malignant from malignant tissue (AUC 0.949, 95% CI 0.891–1.0; ρ = 0.717, P = 6.8 × 10⁻⁹, n = 49). Single-cell localisation: across three blocks of one PVL lesion (GSE196296; 28,543 cells) the between-block effect size was η² = 0.006, excluding linear escalation along the lesion axis, and the signature concentrated in five cell states (Cliff's δ = 0.43–0.89). Drug repositioning: two-direction LINCS reversal returned eight candidates, with PI3K/mTOR inhibitors (torin-1, GSK-1059615, ZSTK-474) prioritised on human oral premalignant lesion precedent (60% histological response to metformin; pS6 reduction correlated with response, P = 0.04), not cell-line reasoning.

We have tried to make our results harder to over-read. The internal cross-validated AUC of 1.0 is labelled an optimistic upper bound, as gene selection and evaluation were not nested. The single-cell η² = 0.006 is reported as an honest null: Kruskal–Wallis returns P = 1.13 × 10⁻¹¹ across 28,543 cells, and we decline to read that as a gradient. GSE85195 is a cross-sectional state-discrimination test in a gingivobuccal complex population, not a prospective transformation prediction. Signature-level enrichment is sparse (one Hallmark term, no transcription-factor or KEGG terms), so mTORC1, PI3K-AKT-mTOR and EMT remain bulk cancer-versus-control findings, not signature mechanisms. Because our signature shares its discovery cohort with Saintigny's, we report zero symbol-level overlap while stating that the cohorts are not independent; independence rests on GSE227919 and GSE85195 alone.

The authors' primary clinical affiliations are in ophthalmology, and this work is a secondary analysis of public transcriptomic datasets rather than a study of our own patients. The study reflects a methodological interest in cross-platform signature transfer and makes no claim of clinical expertise in oral medicine. Oral potentially malignant disorders were chosen because the public data landscape offers an unusually complete pre-malignant-to-malignant continuum — a discovery cohort with transformation outcome, an ordinal-severity RNA-seq cohort, a single lesion sampled along its pathological continuum, and an independent microarray cohort — making it an unusually stringent test bed for a platform-agnostic scoring framework.

The manuscript is original, has not been published in whole or in part, and is not under consideration elsewhere. All authors have read and approved the submission; there is no authorship dispute. No competing interests are declared; this research received no specific funding. Input datasets are public (GSE26549, GSE227919, GSE196296, GSE85195); the signature, per-sample scores and analysis code are archived under a concept DOI (10.5281/zenodo.21823583; version 2, 10.5281/zenodo.23453360) and at https://github.com/PaperValidation/stomatology-pvl. De-identified public data required no ethics approval. Generative AI assistance is declared in the manuscript per journal policy.

We hope the work suits *Oral Oncology*'s readership: it turns an extensively analysed but thinly validated literature into a platform-agnostic score with explicit bounds, and yields testable chemoprevention hypotheses rather than clinical recommendations. We would be glad to respond to reviewer comments and thank you for your consideration.

Yours sincerely,

Wang Yuqian
On behalf of all authors
Shaanxi Eye Hospital, Xi'an People's Hospital (Xi'an Fourth Hospital), Affiliated People's Hospital of Northwest University, Xi'an 710004, China
E-mail: wyq@med.nwu.edu.cn
ORCID: 0009-0007-2276-0085
