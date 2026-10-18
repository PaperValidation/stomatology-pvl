#!/usr/bin/env Rscript
source(file.path(getwd(), "_env_paths.R"))  # paths from .env (single source of truth)
## ============================================================
## PVL Phase 9: SECOND-INDEPENDENT-COHORT cross-sectional discrimination
## Cohort : GSE85195 (Agilent 4x44K, Gingivobuccal complex)
##          groups from !Sample_title: OSCC (tumor) / OPL (leukoplakia) / Normal (n=1)
## Method : IDENTICAL direction-weighted standardized signature score as PVL_04
##          score = safe_scale(expr) %*% w,  w = GSE26549 progression-DEG logFC (FIXED,
##          discovery weights, NO re-estimation on this cohort, NO platform-scale transfer)
## Honest framing: OPL-vs-OSCC is cross-sectional STATE discrimination of
##          premalignant vs malignant, NOT prospective prediction of who transforms.
## ============================================================
suppressMessages({
  library(data.table); library(dplyr); library(pROC); library(ggplot2)
})

RAW  <- file.path(DATASET_PATH, "raw")
FIG  <- file.path(RES, "manuscript_figs")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
cat("=== PVL Phase 9 start ===\n")

## ---------- helpers (same as PVL_04) ----------
safe_scale <- function(mat){            # genes x samples -> samples x genes z-scored per gene
  m <- t(mat)
  sds <- apply(m, 2, sd, na.rm = TRUE)
  if (any(sds > 1e-8)) m[, sds > 1e-8] <- scale(m[, sds > 1e-8, drop = FALSE])
  m[, sds <= 1e-8] <- 0
  m
}
parse_gse_matrix <- function(txt){
  lines <- readLines(txt)
  beg <- which(grepl("!series_matrix_table_begin", lines))[1] + 1
  end <- which(grepl("!series_matrix_table_end",   lines))[1]
  block <- lines[beg:end]
  hdr <- gsub('"', '', strsplit(block[1], "\t")[[1]])
  m <- fread(text = block[-1], sep = "\t", header = FALSE, fill = TRUE)
  rn <- as.character(m[[1]]); m <- as.matrix(m[, -1, with = FALSE])
  rownames(m) <- rn; colnames(m) <- hdr[-1]; m
}

## ---------- 1. discovery signature weights w (FIXED to 101-gene signature) ----------
deg <- fread(file.path(RES, "GSE26549_progression_DEG.csv"))
sig <- unique(deg[adj.P.Val < 0.05 & abs(logFC) > 0.5]$gene)
sig <- sig[!is.na(sig)]
## FIX: restrict weights to the discovery signature ONLY (matches PVL_04 headline;
## the DEG csv contains the full contrast gene list, not just the 101-gene signature)
w_sig <- setNames(deg$logFC, deg$gene)
w_sig <- w_sig[!is.na(w_sig)]
w_sig <- w_sig[names(w_sig) %in% sig]
cat("discovery signature genes (adj.P<0.05 & |logFC|>0.5):", length(sig),
    " ; fixed signature weights w:", length(w_sig), "\n")

## ---------- 2. GSE85195 matrix + sample titles ----------
mat_file <- file.path(RAW, "GSE85195", "GSE85195_series_matrix.txt")
L <- readLines(mat_file)
title_line <- L[grep("^!Sample_title", L)[1]]
titles <- gsub('"', '', strsplit(title_line, "\t")[[1]][-1])
expr <- parse_gse_matrix(mat_file)                 # probes x samples (numeric)
cat("GSE85195 matrix:", nrow(expr), "probes x", ncol(expr), "samples\n")
stopifnot(length(titles) == ncol(expr))

## ---------- 3. probe -> gene (Agilent GPL6480) ----------
pg <- fread(file.path(RAW, "GPL6480_probe2gene.tsv"))
pg <- pg[!is.na(gene) & gene != ""]
mt <- expr[rownames(expr) %in% pg$probe, , drop = FALSE]
pgm <- pg[match(rownames(mt), pg$probe), ]
genes <- unique(pgm$gene)
exprG <- sapply(genes, function(g){
  idx <- which(pgm$gene == g)
  if (length(idx) == 1) mt[idx, ] else apply(mt[idx, , drop = FALSE], 2, max)
})
exprG <- as.matrix(exprG); exprG <- t(exprG)        # genes x samples
cat("collapsed", nrow(mt), "probes ->", length(genes), "genes\n")

## ---------- 4. intersect signature with this cohort, fix weights ----------
genes_avail <- intersect(names(w_sig), rownames(exprG))
w2 <- w_sig[genes_avail]
cat("signature genes present in GSE85195:", length(w2),
    " (mapping rate ", round(length(w2)/length(w_sig)*100, 1), "%)\n")

## ---------- 5. classify groups from sample titles ----------
group <- ifelse(grepl("OSCC", titles), "OSCC",
         ifelse(grepl("OPL",  titles), "OPL",
         ifelse(grepl("Control", titles), "Normal", "Other")))
names(group) <- colnames(expr)
tab <- table(group)
cat("group sizes:", paste(names(tab), tab, collapse = " / "), "\n")

## ---------- 6. direction-weighted signature score ----------
Z    <- safe_scale(exprG[names(w2), , drop = FALSE])   # samples x genes
score <- as.numeric(Z %*% w2)
names(score) <- colnames(exprG)

## ---------- 7. OPL vs OSCC discrimination (primary, cross-sectional) ----------
idx_disc <- group %in% c("OPL", "OSCC")
yd <- as.integer(group[idx_disc] == "OSCC")        # OSCC = positive
sd_ <- score[idx_disc]
roc_d <- roc(yd, sd_)
auc_d <- as.numeric(auc(roc_d))
ci_d  <- ci(roc_d)
cat("OPL-vs-OSCC discrimination AUC =", round(auc_d, 3),
    " 95%CI", round(min(ci_d),3), "-", round(max(ci_d),3), "\n")

## ---------- 8. severity ordinal (OPL < OSCC) corroboration ----------
ord <- ifelse(group == "OPL", 0, ifelse(group == "OSCC", 1, NA))
keep <- !is.na(ord)
rho_sev <- as.numeric(cor.test(score[keep], ord[keep], method = "spearman")$estimate)
p_sev   <- as.numeric(cor.test(score[keep], ord[keep], method = "spearman")$p.value)
cat("Severity ordinal (OPL<OSCC) Spearman rho =", round(rho_sev, 3), " p =", signif(p_sev, 2), "\n")

## mean score per group
ms <- tapply(score, group, mean)
cat("mean score: ", paste(names(ms), round(ms,3), collapse = " / "), "\n")

## ---------- 9. outputs ----------
out_df <- data.frame(sample = colnames(expr), gsm = colnames(expr),
                     group = group, title = titles, score = score,
                     stringsAsFactors = FALSE)
write.csv(out_df, file.path(FIG, "PVL9_GSE85195_scores.csv"), row.names = FALSE)
write.csv(data.frame(metric = c("n_OPL","n_OSCC","n_Normal",
                                "sig_genes_total","sig_genes_mapped","mapping_rate_pct",
                                "OPL_vs_OSCC_AUC","AUC_CI_low","AUC_CI_high",
                                "severity_rho_OPL_OSCC","severity_p",
                                "mean_score_OPL","mean_score_OSCC","mean_score_Normal"),
                     value = c(sum(group=="OPL"), sum(group=="OSCC"), sum(group=="Normal"),
                               length(sig), length(w2), round(length(w2)/length(sig)*100,1),
                               auc_d, min(ci_d), max(ci_d), rho_sev, p_sev,
                               ms["OPL"], ms["OSCC"], ms["Normal"])),
          file.path(FIG, "PVL9_summary.csv"), row.names = FALSE)

## ROC plot
png(file.path(FIG, "PVL9_ROC_OPL_vs_OSCC.png"), width = 700, height = 650, res = 120)
plot(roc_d, main = "GSE85195: OPL vs OSCC\ndirection-weighted signature score",
     col = "#1f78b4", lwd = 2)
abline(0, 1, lty = 2, col = "grey")
legend("bottomright", sprintf("AUC = %.3f", auc_d), bty = "n")
dev.off()

## boxplot by group
bp <- data.frame(group = factor(group, c("Normal","OPL","OSCC")), score = score)
png(file.path(FIG, "PVL9_boxplot_score.png"), width = 700, height = 600, res = 120)
ggplot(bp, aes(x = group, y = score, fill = group)) +
  geom_boxplot(alpha = .7) + geom_jitter(width = .15, size = 1.5) +
  scale_fill_manual(values = c("Normal"="#999999","OPL"="#1f78b4","OSCC"="#e31a1c")) +
  labs(title = "Direction-weighted signature score by group (GSE85195)",
       y = "signature score", x = "") + theme_minimal()
dev.off()

cat("\n=== PVL Phase 9 DONE ===\n")
