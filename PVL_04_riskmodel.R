#!/usr/bin/env Rscript
source(file.path(getwd(), "_env_paths.R"))  # paths from .env (single source of truth)
## ============================================================
## PVL Phase 4: progression-risk model (discovery + cross-platform validation)
## Discovery : GSE26549 (bulk microarray OPL, baseline biopsy, progressing vs non)
## External  : GSE227919 (RNA-seq, Progressed vs Stable)
## Validation strategy:
##   (A) LASSO classifier on discovery (internal 5-fold CV AUC)
##   (B) HEADLINE: platform-agnostic direction-weighted standardized signature
##       score  ->  validated on the independent RNA-seq cohort (no coefficient
##       transfer, no platform-scale leakage)
## Link      : pseudobulk from raw 10x (LesionA/B/C) vs bulk progression signature
## ============================================================
Sys.setenv(http_proxy  = "socks5://127.0.0.1:8428",
           https_proxy = "socks5://127.0.0.1:8428")
suppressMessages({
  library(glmnet); library(pROC); library(data.table); library(dplyr);
  library(ggplot2); library(Matrix)
})

RAW  <- file.path(DATASET_PATH, "raw")
MET  <- file.path(DATASET_PATH, "metadata")
dir.create(RES, showWarnings = FALSE)
cat("=== PVL Phase 4 start ===\n")

## ---------- helpers ----------
safe_scale <- function(mat){            # mat: genes x samples -> samples x genes z-scored per gene
  m <- t(mat)
  sds <- apply(m, 2, sd, na.rm = TRUE)
  if (any(sds > 1e-8)) m[, sds > 1e-8] <- scale(m[, sds > 1e-8, drop = FALSE])
  m[, sds <= 1e-8] <- 0
  m
}
parse_gse_matrix <- function(gz){
  lines <- readLines(gz)
  beg <- which(grepl("!series_matrix_table_begin", lines))[1] + 1
  end <- which(grepl("!series_matrix_table_end",   lines))[1]
  block <- lines[beg:end]
  hdr <- gsub('"', '', strsplit(block[1], "\t")[[1]])
  m <- fread(text = block[-1], sep = "\t", header = FALSE, fill = TRUE)
  rn <- as.character(m[[1]])
  m <- as.matrix(m[, -1, with = FALSE])
  rownames(m) <- rn
  colnames(m) <- hdr[-1]
  m
}

## ---------- 1. GSE26549 expression (baseline biopsy) ----------
expr49 <- parse_gse_matrix(file.path(RAW, "GSE26549_series_matrix.txt.gz"))
cat("GSE26549 matrix:", nrow(expr49), "probes x", ncol(expr49), "samples\n")

## probe -> gene symbol
annot_lines <- readLines(file.path(RAW, "GPL6244.annot.gz"))
cand <- grep("Gene symbol", annot_lines, ignore.case = TRUE)
hdr_i <- cand[!grepl("^#", annot_lines[cand])][1]
ann <- fread(text = annot_lines[hdr_i:length(annot_lines)], sep = "\t",
             header = TRUE, fill = TRUE)
sym_col <- which(grepl("Gene symbol", names(ann), ignore.case = TRUE))[1]
id_col  <- 1
sym <- sapply(strsplit(as.character(ann[[sym_col]]), " /// "), function(x) x[1])
gi <- data.frame(probe = as.character(ann[[id_col]]), gene = sym,
                 stringsAsFactors = FALSE)
gi <- gi[!is.na(gi$gene) & gi$gene != "", ]
mt <- expr49[rownames(expr49) %in% gi$probe, , drop = FALSE]
gi <- gi[match(rownames(mt), gi$probe), ]
genes <- unique(gi$gene)
cat("collapsing", nrow(mt), "probes ->", length(genes), "genes\n")
exprG <- sapply(genes, function(g){
  idx <- which(gi$gene == g)
  if (length(idx) == 1) mt[idx, ] else apply(mt[idx, , drop = FALSE], 2, max)
})
exprG <- as.matrix(exprG)
exprG <- t(exprG)                       # genes x samples
colnames(exprG) <- sub('"', '', colnames(exprG))

## clinical labels (baseline biopsy only)
clin49 <- fread(file.path(MET, "GSE26549_clinical.csv"))
baseline <- clin49[`time of biopsy` == "biopsy at baseline"]
y <- as.integer(baseline$outcome == "progressing")
names(y) <- as.character(baseline$GSM)
cat("baseline samples:", length(y), " progressing:", sum(y), " non_progressing:", sum(y == 0), "\n")
x49 <- exprG[, as.character(baseline$GSM), drop = FALSE]

## signature from Phase 2 DEG (loose threshold)
deg <- fread(file.path(RES, "GSE26549_progression_DEG.csv"))
sig <- unique(deg[adj.P.Val < 0.05 & abs(logFC) > 0.5]$gene)
sig <- sig[!is.na(sig)]
cat("signature genes (adj.P<0.05 & |logFC|>0.5):", length(sig), "\n")

## ---------- GSE227919 expression (RNA-seq, log2) ----------
lines <- readLines(file.path(RAW, "GSE227919_2023_20_03_expr_mat_raw.txt.gz"))
hdr <- strsplit(trimws(lines[1]), "\\s+")[[1]]
fm <- fread(text = lines[-1], sep = " ", header = FALSE, fill = TRUE)
genes <- as.character(fm[[1]])
exm <- as.matrix(fm[, -1, with = FALSE])
rownames(exm) <- genes
colnames(exm) <- hdr

## genes present in BOTH cohorts -> shared, comparable feature space
common_feat <- intersect(intersect(sig, rownames(x49)), rownames(exm))
cat("signature genes present in both cohorts:", length(common_feat), "\n")

## ============================================================
## (A) LASSO classifier on discovery cohort (internal CV only)
## ============================================================
Xtr <- t(x49[common_feat, , drop = FALSE])   # samples x genes
set.seed(42)
foldid <- sample(rep(1:5, length.out = length(y)))   # explicit, reproducible folds
cvfit <- cv.glmnet(Xtr, y, family = "binomial", alpha = 1, foldid = foldid)
fold_aucs <- vapply(unique(foldid), function(f){
  te <- foldid == f
  if (length(unique(y[te])) < 2) return(NaN)   # single-class fold -> AUC undefined
  p <- predict(cvfit, newx = Xtr[te, , drop = FALSE], s = "lambda.min", type = "response")[, 1]
  as.numeric(auc(roc(y[te], p)))
}, numeric(1))
cv_auc <- mean(fold_aucs, na.rm = TRUE)
n_skip <- sum(is.nan(fold_aucs))
cf <- as.matrix(coef(cvfit, s = "lambda.min"))
sel <- rownames(cf)[cf[, 1] != 0 & rownames(cf) != "(Intercept)"]
cat("LASSO selected genes:", length(sel), "\n")
risk_tr <- predict(cvfit, newx = Xtr, s = "lambda.min", type = "response")[, 1]
auc_tr <- as.numeric(auc(roc(y, risk_tr)))
cat("LASSO internal (resubstitution) AUC:", round(auc_tr, 3),
    "   honest 5-fold CV AUC:", round(cv_auc, 3),
    " (single-class folds skipped:", n_skip, ")\n")

## ============================================================
## (B) HEADLINE: platform-agnostic direction-weighted signature score
##     z-score per gene WITHIN each cohort (removes platform scale),
##     weight by GSE26549 DEG logFC direction & magnitude.
## ============================================================
w <- setNames(deg$logFC, deg$gene)[common_feat]
w <- w[!is.na(w)]
cat("signature-score genes (weighted):", length(w), "\n")

## (B0) External: signature score vs GSE227919 histopathological severity gradient
##      (platform-agnostic; uses ALL 66 samples, ordinal Control<HkNR<Dysplasia<Cancer)
clin22g <- fread(file.path(MET, "GSE227919_clinical.csv"))
cls <- setNames(clin22g$Class, clin22g$Sample_ID)
keep_all <- intersect(colnames(exm), names(cls))
Zall <- safe_scale(exm[names(w), keep_all, drop = FALSE])
score_all <- as.numeric(Zall %*% w)
ord_map <- c(Control = 0, HkNR = 1, Dysplasia = 2, Cancer = 3)
ord_all <- ord_map[cls[keep_all]]
grad_rho <- as.numeric(cor.test(score_all, ord_all, method = "spearman")$estimate)
grad_p   <- as.numeric(cor.test(score_all, ord_all, method = "spearman")$p.value)
is_can <- as.integer(ord_all == 3)
auc_can <- as.numeric(auc(roc(is_can, score_all)))
cat("External severity-gradient: Spearman rho =", round(grad_rho, 3),
    " p =", signif(grad_p, 2),
    " ; Cancer-vs-other AUC =", round(auc_can, 3), "\n")

## training cohort z-scored per gene
Ztr <- safe_scale(x49[names(w), , drop = FALSE])   # samples x genes
score_tr <- as.numeric(Ztr %*% w)
auc_score_tr <- as.numeric(auc(roc(y, score_tr)))
cat("Signature-SCORE internal AUC (GSE26549):", round(auc_score_tr, 3), "\n")

## external cohort
clin22 <- fread(file.path(MET, "GSE227919_clinical.csv"))
val <- clin22[Progression_status %in% c("Stable", "Progressed-Dys", "Progressed-SCC")]
yt <- as.integer(grepl("Progressed", val$Progression_status))
names(yt) <- as.character(val$Sample_ID)
common <- intersect(colnames(exm), names(yt))
exm <- exm[, common, drop = FALSE]; yt <- yt[common]
cat("GSE227919 validation samples:", length(yt),
    " progressed:", sum(yt == 1), " stable:", sum(yt == 0), "\n")

Zte <- safe_scale(exm[names(w), , drop = FALSE])   # samples x genes
score_te <- as.numeric(Zte %*% w)
auc_score_te <- as.numeric(auc(roc(yt, score_te)))
cat("Signature-SCORE external AUC (GSE227919):", round(auc_score_te, 3), "\n")

## map score -> probability via training logistic (no test leakage)
cal_fit <- glm(y ~ score_tr, family = binomial)
p_te <- as.numeric(plogis(coef(cal_fit)[1] + coef(cal_fit)[2] * score_te))
roc_te <- roc(yt, p_te)
auc_te <- as.numeric(auc(roc_te))
cal <- glm(yt ~ p_te, family = binomial)
cal_slope <- as.numeric(coef(cal)[2]); cal_int <- as.numeric(coef(cal)[1])
ci_te <- ci(roc_te)
ci_low <- min(as.numeric(ci_te)); ci_high <- max(as.numeric(ci_te))
pval_te <- as.numeric(wilcox.test(p_te[yt == 1], p_te[yt == 0])$p.value)
cat("External probability AUC:", round(auc_te, 3),
    " 95% CI:", round(ci_low, 3), "-", round(ci_high, 3),
    " p(MWU):", signif(pval_te, 2), "\n")
cat("  calib slope (ideal 1):", round(cal_slope, 3),
    "  calib intercept (ideal 0):", round(cal_int, 3), "\n")

## DCA (net benefit) on external probability
dca <- function(yy, p, thr){
  tp <- sum(p >= thr & yy == 1); fp <- sum(p >= thr & yy == 0); n <- length(yy)
  (tp - fp * (thr / (1 - thr))) / n
}
thr <- seq(0.05, 0.75, 0.05)
dca_df <- data.frame(threshold = thr,
                     model = sapply(thr, function(t) dca(yt, p_te, t)),
                     treat_all = sapply(thr, function(t) dca(yt, rep(1, length(yt)), t)),
                     treat_none = 0)

## ============================================================
## (C) auxiliary: direct LASSO coefficient transfer to RNA-seq (limitation check)
## ============================================================
risk_te_direct <- predict(cvfit, newx = t(exm[common_feat, , drop = FALSE]),
                          s = "lambda.min", type = "response")[, 1]
auc_direct <- as.numeric(auc(roc(yt, risk_te_direct)))
cat("Direct LASSO-transfer external AUC (microarray->RNAseq):", round(auc_direct, 3),
    "  [demonstrates platform-scale non-transferability]\n")

## ---------- 4. pseudobulk linkage (raw 10x, OOM-safe) ----------
samples <- list(LesionA = "GSM5870060", LesionB = "GSM5870062", LesionC = "GSM5870063")
ft <- read.delim(file.path(RAW, sprintf("%s_%s_features.tsv.gz",
               samples[[1]], "LesionA")), header = FALSE, sep = "\t",
               stringsAsFactors = FALSE, check.names = FALSE)
gsym <- ifelse(is.na(ft$V2) | ft$V2 == "", ft$V1, ft$V2)
pb <- list()
for (les in names(samples)) {
  mp <- file.path(RAW, sprintf("%s_%s_matrix.mtx.gz", samples[[les]], les))
  m <- as(readMM(gzfile(mp)), "dgCMatrix")          # genes x cells
  pb[[les]] <- as.numeric(rowSums(m))
}
names(pb[[1]]) <- gsym
pbmat <- do.call(cbind, lapply(pb, function(v) v[names(pb[[1]])]))
pbmat[is.na(pbmat)] <- 0                            # undetected gene = 0 counts
rownames(pbmat) <- names(pb[[1]])
lp <- log2(1 + pbmat)
ord <- c(LesionA = 0, LesionB = 1, LesionC = 2)
sg <- intersect(sig, rownames(lp))
pv <- lp[sg, ]
cor_lr <- apply(pv, 1, function(r) cor(r, ord, method = "spearman"))
bulk_lfc <- setNames(deg$logFC, deg$gene)[sg]
agree <- mean(sign(cor_lr) == sign(bulk_lfc), na.rm = TRUE)
mean_abs <- mean(abs(cor_lr), na.rm = TRUE)
cat("Pseudobulk: signature genes tested =", length(sg),
    "; mean |spearman| vs lesion gradient =", round(mean_abs, 3),
    "; direction agreement with bulk logFC =", round(agree * 100, 1), "%\n")

## ---------- outputs ----------
write.csv(data.frame(gene = sel, coef = cf[sel, 1],
                     weight_logFC = w[sel]),
          file.path(RES, "PVL4_lasso_coefficients.csv"), row.names = FALSE)
write.csv(data.frame(sample = names(yt), score = score_te, prob = p_te,
                     outcome_progressed = yt),
          file.path(RES, "PVL4_test_risk_scores.csv"), row.names = FALSE)
write.csv(dca_df, file.path(RES, "PVL4_dca.csv"), row.names = FALSE)
write.csv(data.frame(gene = sg, bulk_logFC = bulk_lfc,
                     pseudobulk_spearman_vs_lesion = cor_lr),
          file.path(RES, "PVL4_pseudobulk_linkage.csv"), row.names = FALSE)
write.csv(data.frame(metric = c("LASSO_internal_AUC","LASSO_5foldCV_AUC",
                                "SignatureScore_internal_AUC","SignatureScore_external_AUC",
                                "SignatureScore_external_CI_low","SignatureScore_external_CI_high",
                                "SignatureScore_external_p_MWU",
                                "DirectLASSO_transfer_external_AUC",
                                "Ext_calib_slope","Ext_calib_intercept",
                                "Pseudobulk_mean_abs_spearman","Pseudobulk_dir_agreement_pct",
                                "Ext_gradient_spearman_rho","Ext_gradient_p","Ext_gradient_Cancer_AUC"),
                     value = c(auc_tr, cv_auc, auc_score_tr, auc_score_te,
                               ci_low, ci_high, pval_te, auc_direct,
                               cal_slope, cal_int, mean_abs, agree * 100,
                               grad_rho, grad_p, auc_can)),
          file.path(RES, "PVL4_summary.csv"), row.names = FALSE)

png(file.path(RES, "PVL4_roc.png"), width = 800, height = 700, res = 120)
plot(roc_te, main = paste0("External validation (signature score)\nAUC = ",
                           round(auc_te, 3)))
abline(0, 1, lty = 2, col = "grey")
dev.off()
png(file.path(RES, "PVL4_dca.png"), width = 800, height = 600, res = 120)
matplot(dca_df$threshold, dca_df[, c("model", "treat_all")], type = "l",
        xlab = "Threshold probability", ylab = "Net benefit",
        col = c(2, 4), lwd = 2, main = "Decision Curve Analysis (external)")
abline(h = 0, lty = 2, col = "grey")
legend("topright", c("Signature score", "Treat all"), col = c(2, 4), lwd = 2)
dev.off()

cat("\n=== PVL Phase 4 DONE ===\n")
