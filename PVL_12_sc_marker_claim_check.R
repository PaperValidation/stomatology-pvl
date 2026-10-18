#!/usr/bin/env Rscript
# PVL_12: marker-defined check of the Treg / plasma / mast claim (claim §3.3)
# READ-ONLY. Only reads a small gene panel from the counts layer.
suppressPackageStartupMessages(library(Seurat))
source(file.path(getwd(), "_env_paths.R"))

rds_path <- file.path(RES, "PVL_sc_combined.rds")
out_fig  <- file.path(getwd(), "results", "manuscript_figs")
out_data <- RES

cat("=== load ===\n"); t0 <- Sys.time()
obj <- readRDS(rds_path)
cat("loaded in", round(as.numeric(difftime(Sys.time(), t0, units="secs")),1), "s\n")
md <- obj@meta.data
md$lesion <- as.character(md$lesion)
md$celltype <- as.character(md$celltype)

panel <- c("PTPRC","CD3D","CD3E","CD2","TRAC",          # pan-T
           "FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18",  # Treg
           "MZB1","JCHAIN","IGHG1","IGKC","SDC1","XBP1","DERL3", # plasma
           "TPSAB1","TPSB2","CPA3","KIT","MS4A2","HDC",  # mast
           "KRT10","KRT14","KRT1","KRT16","SBSN","LOR",  # keratinocyte
           "CD19","MS4A1","CD79A")                       # B
genes <- intersect(panel, rownames(obj))
cat("missing from object:", paste(setdiff(panel, rownames(obj)), collapse=", "), "\n")

cnt <- LayerData(obj, assay = "RNA", layer = "counts")
m <- as(cnt[genes, , drop = FALSE], "dgCMatrix")
rm(cnt); rm(obj); invisible(gc())
cat("panel matrix:", paste(dim(m), collapse=" x "), "\n")

expr <- m > 0   # detected / not
les  <- md$lesion
lev  <- c("LesionA","LesionB","LesionC")
N    <- sapply(lev, function(l) sum(les == l))

det_rate <- function(g, subset_idx = rep(TRUE, ncol(m))) {
  sapply(lev, function(l) {
    idx <- subset_idx & les == l
    if (sum(idx) == 0) return(NA_real_)
    100 * mean(expr[g, idx])
  })
}
mean_pos <- function(g, subset_idx = rep(TRUE, ncol(m))) {
  sapply(lev, function(l) {
    idx <- subset_idx & les == l
    if (sum(idx) == 0) return(NA_real_)
    mean(m[g, idx])
  })
}

cat("\n=== detection rate (%) per lesion, all cells ===\n")
rows <- lapply(genes, function(g) data.frame(
  gene = g,
  pct_LesionA = round(det_rate(g)["LesionA"], 3),
  pct_LesionB = round(det_rate(g)["LesionB"], 3),
  pct_LesionC = round(det_rate(g)["LesionC"], 3),
  mean_A = round(mean_pos(g)["LesionA"], 4),
  mean_B = round(mean_pos(g)["LesionB"], 4),
  mean_C = round(mean_pos(g)["LesionC"], 4),
  stringsAsFactors = FALSE))
dr <- do.call(rbind, rows)
print(dr, row.names = FALSE)
write.csv(dr, file.path(out_fig, "SC_marker_detection_by_lesion.csv"), row.names = FALSE)

# ---- marker-defined states ----
E <- function(g) expr[g, ]
V <- function(g) as.numeric(m[g, ])          # raw counts
isT   <- E("PTPRC") & E("CD3D") & E("CD3E")
isTreg <- isT & E("FOXP3") & (E("IL2RA") | E("CTLA4") | E("TNFRSF18"))
isKera   <- E("KRT10") | E("KRT14")
# NOTE: detection>0 is strongly inflated by ambient/soup RNA (KRT14, KRT16, IGKC
# are detected in ~100% of cells), so state definitions below use count thresholds.
isTreg_strict <- E("CD3D") & E("CD3E") & (V("FOXP3") >= 2)
isPlasma1 <- E("MZB1") & E("JCHAIN")                       # >0
isPlasma2 <- (V("MZB1") >= 2) & E("JCHAIN")                # MZB1>=2
isPlasma3 <- (V("MZB1") >= 2) & (V("JCHAIN") >= 2)         # both >=2
isPlasma4 <- (V("MZB1") >= 3) & (V("JCHAIN") >= 3)         # strict
isPlasma5 <- (V("MZB1") >= 2) & (V("JCHAIN") >= 2) & !isT  # strict + not T
isMast1   <- E("TPSAB1") & E("TPSB2") & E("CPA3")          # >0
isMast2   <- (V("TPSAB1") >= 2) & (V("TPSB2") >= 1)        # thresholded
isMast3   <- (V("TPSAB1") >= 2) & (V("TPSB2") >= 1) & E("CPA3") & E("KIT")

states <- list(T_cells_marker = isT, Treg_marker = isTreg, Treg_marker_strict = isTreg_strict,
               plasma_marker_gt0 = isPlasma1, plasma_MZB1ge2 = isPlasma2,
               plasma_both_ge2 = isPlasma3, plasma_both_ge3 = isPlasma4,
               plasma_both_ge2_notT = isPlasma5,
               mast_gt0 = isMast1, mast_TPSAB1ge2 = isMast2, mast_strict = isMast3)
cat("\n=== marker-defined state counts / % within lesion block ===\n")
out <- list()
for (nm in names(states)) {
  v <- states[[nm]]
  n <- sapply(lev, function(l) sum(v & les == l))
  p <- 100 * n / N
  # 2x3 test: state vs rest across lesions
  mat <- rbind(n, N - n)
  ch  <- suppressWarnings(chisq.test(mat, correct = FALSE))
  pv  <- tryCatch(suppressWarnings(prop.test(x = as.numeric(n), n = as.numeric(N), correct = FALSE))$p.value,
                  error = function(e) NA_real_)
  chi2 <- as.numeric(ch$statistic); tot <- sum(mat)
  V <- sqrt(chi2 / (tot * (min(dim(mat)) - 1)))
  out[[nm]] <- data.frame(state = nm,
                          n_A = n[1], n_B = n[2], n_C = n[3],
                          pct_A = round(p[1],3), pct_B = round(p[2],3), pct_C = round(p[3],3),
                          max_minus_min = round(max(p) - min(p), 3),
                          fold_C_over_max_AB = round(p[3]/max(p[1],p[2],1e-9), 2),
                          p_value = signif(pv, 4), cramers_V = round(V, 4),
                          stringsAsFactors = FALSE)
}
sm <- do.call(rbind, out); sm$FDR <- signif(p.adjust(sm$p_value, "BH"), 4)
print(sm, row.names = FALSE)
write.csv(sm, file.path(out_fig, "SC_marker_state_by_lesion_stats.csv"), row.names = FALSE)

# ---- denominator-restricted: Treg as % of marker-defined T cells ----
cat("\n=== Treg as fraction of marker-defined T cells (denominator = CD3D+CD3E+PTPRC+) ===\n")
nT <- sapply(lev, function(l) sum(isT & les == l))
nTr <- sapply(lev, function(l) sum(isTreg & les == l))
cat(sprintf("  T cells: A=%d B=%d C=%d\n", nT[1], nT[2], nT[3]))
cat(sprintf("  Treg/T:   A=%.2f%% B=%.2f%% C=%.2f%%\n",
            100*nTr[1]/nT[1], 100*nTr[2]/nT[2], 100*nTr[3]/nT[3]))
pt2 <- suppressWarnings(prop.test(x = as.numeric(nTr), n = as.numeric(nT), correct = FALSE))
cat("  prop.test p =", signif(pt2$p.value, 4), "\n")

# ---- within-label check: cluster 6 (FOXP3-high) occupancy ----
cat("\n=== cluster 6 / 19 / 10 / 18 occupancy by lesion (from meta) ===\n")
if ("seurat_clusters" %in% colnames(md)) {
  cl <- as.character(md$seurat_clusters)
  for (k in c("6","19","10","18")) {
    tb <- table(factor(les[cl == k], levels = lev))
    cat(sprintf("  cluster %-3s n=%4d  A=%4d B=%4d C=%4d   (%%of cluster in C = %.1f%%)\n",
                k, sum(tb), tb[1], tb[2], tb[3], 100*tb[3]/sum(tb)))
  }
}

cat("\n=== DONE ===\n")
