#!/usr/bin/env Rscript
# PVL_11: celltype x lesion cross-tabs + stats for the manuscript claim audit
# READ-ONLY on the Seurat rds. Adds new CSVs only.
suppressPackageStartupMessages(library(Seurat))
source(file.path(getwd(), "_env_paths.R"))

rds_path <- file.path(RES, "PVL_sc_combined.rds")
out_fig  <- file.path(getwd(), "results", "manuscript_figs")
out_data <- RES
dir.create(out_fig,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)

cat("=== loading rds (read-only) ===\n")
t0 <- Sys.time()
obj <- readRDS(rds_path)
cat("loaded in", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
cat("class:", class(obj), "\n")
cat("dim:", paste(dim(obj), collapse = " x "), "\n")

md <- obj@meta.data
cat("meta.data dim:", paste(dim(md), collapse = " x "), "\n")
cat("meta cols:", paste(colnames(md), collapse = ", "), "\n")
print(utils::head(md, 3))

# free the big object asap
rm(obj); invisible(gc())

stopifnot(all(c("celltype", "lesion") %in% colnames(md)))

md$celltype <- as.character(md$celltype)
md$lesion   <- as.character(md$lesion)

cat("\n=== table(celltype) full ===\n")
print(table(md$celltype, useNA = "ifany"))
cat("\n=== table(lesion) ===\n")
print(table(md$lesion, useNA = "ifany"))
cat("\nNA celltype n =", sum(is.na(md$celltype)), "\n")

# drop NA celltype for cross-tabs, but report
md_c <- md[!is.na(md$celltype), ]

lev_les <- intersect(c("LesionA", "LesionB", "LesionC"), unique(md_c$lesion))
cat("lesion order used:", paste(lev_les, collapse = ", "), "\n")

# ---- T1: celltype x lesion counts ----
ct <- table(factor(md_c$celltype, levels = sort(unique(md_c$celltype))),
            factor(md_c$lesion,   levels = lev_les))
ct_df <- as.data.frame.matrix(ct, stringsAsFactors = FALSE)
ct_df <- cbind(celltype = rownames(ct_df), ct_df, row.names = NULL, stringsAsFactors = FALSE)
ct_df$n_total <- rowSums(ct)
ct_df <- ct_df[order(-ct_df$n_total), ]
write.csv(ct_df, file.path(out_fig, "SC_celltype_by_lesion_counts.csv"), row.names = FALSE)
cat("\n[written] SC_celltype_by_lesion_counts.csv :",
    nrow(ct_df), "x", ncol(ct_df), "\n")
print(ct_df)

# ---- T1: within-lesion proportions (column normalised) ----
props <- sweep(ct, 2, colSums(ct), "/") * 100
prop_df <- as.data.frame.matrix(props, stringsAsFactors = FALSE)
prop_df <- cbind(celltype = rownames(prop_df), prop_df, row.names = NULL, stringsAsFactors = FALSE)
prop_df <- prop_df[match(ct_df$celltype, prop_df$celltype), ]
prop_df[lev_les] <- lapply(prop_df[lev_les], function(x) round(as.numeric(x), 4))
prop_df$n_total <- ct_df$n_total
write.csv(prop_df, file.path(out_fig, "SC_celltype_by_lesion_props.csv"), row.names = FALSE)
cat("[written] SC_celltype_by_lesion_props.csv :",
    nrow(prop_df), "x", ncol(prop_df), "\n")
cat("column sums (should be ~100):", paste(round(colSums(props), 4), collapse = " / "), "\n")
print(prop_df)

# ---- cluster x lesion x celltype detail ----
if ("seurat_clusters" %in% colnames(md)) {
  cl <- md[, c("seurat_clusters", "lesion", "celltype")]
  cl$seurat_clusters <- as.character(cl$seurat_clusters)
  cl$lesion <- as.character(cl$lesion)
  cl$celltype[is.na(cl$celltype)] <- "NA"
  agg <- aggregate(list(n_cells = rep(1, nrow(cl))),
                   by = list(cluster = cl$seurat_clusters,
                             lesion  = cl$lesion,
                             celltype = cl$celltype), FUN = sum)
  numsort <- function(x) as.integer(as.character(x))
  agg$cluster <- as.character(agg$cluster)
  agg <- agg[order(numsort(agg$cluster), agg$lesion), ]
  write.csv(agg, file.path(out_data, "SC_cluster_lesion_celltype_counts.csv"), row.names = FALSE)
  cat("[written] SC_cluster_lesion_celltype_counts.csv :",
      nrow(agg), "x", ncol(agg), "\n")
  # also a wide cluster x lesion table for readability (new file, not overwriting)
  cl_w <- as.data.frame.matrix(table(cl$seurat_clusters, cl$lesion))
  cl_w$cluster <- rownames(cl_w)
  cl_w <- cl_w[, c("cluster", lev_les)]
  cl_w$celltype <- aggregate(cl$celltype, list(cl$seurat_clusters), function(z) z[1])$x[match(cl_w$cluster, aggregate(cl$celltype, list(cl$seurat_clusters), function(z) z[1])$Group.1)]
  cl_w <- cl_w[order(numsort(cl_w$cluster)), ]
  write.csv(cl_w, file.path(out_data, "SC_cluster_lesion_wide.csv"), row.names = FALSE)
  cat("[written] SC_cluster_lesion_wide.csv :", nrow(cl_w), "x", ncol(cl_w), "\n")
  print(cl_w)
} else {
  cat("!! no seurat_clusters column\n")
}

# ---- T2: per-celltype test across lesions + Cramer's V ----
N_les <- colSums(ct)
res <- list()
for (i in seq_len(nrow(ct))) {
  ctv <- rownames(ct)[i]
  obs <- as.numeric(ct[i, ])                       # this celltype per lesion
  tot <- as.numeric(N_les)
  mat <- rbind(obs, tot - obs)                     # 2 x 3
  n   <- sum(mat)
  pct <- 100 * obs / tot

  ch <- suppressWarnings(tryCatch(chisq.test(mat, correct = FALSE),
                                  error = function(e) NULL))
  exp_min <- if (!is.null(ch)) min(ch$expected) else NA_real_
  if (!is.null(ch) && !is.na(exp_min) && exp_min >= 5) {
    pt <- suppressWarnings(prop.test(x = obs, n = tot, correct = FALSE))
    p <- pt$p.value; test <- "prop.test(2x3, no continuity)"
  } else {
    ft <- tryCatch(fisher.test(mat, simulate.p.value = TRUE, B = 20000),
                   error = function(e) NULL)
    if (!is.null(ft)) { p <- ft$p.value; test <- "fisher.test(sim, B=20000)" }
    else { p <- NA_real_; test <- "failed" }
  }
  chi2 <- if (!is.null(ch)) as.numeric(ch$statistic) else NA_real_
  k <- min(nrow(mat), ncol(mat))
  V <- if (!is.na(chi2) && n > 0 && k > 1) sqrt(chi2 / (n * (k - 1))) else NA_real_

  res[[i]] <- data.frame(
    celltype   = ctv,
    n_total    = sum(obs),
    pct_LesionA = round(pct[match("LesionA", lev_les)], 4),
    pct_LesionB = round(pct[match("LesionB", lev_les)], 4),
    pct_LesionC = round(pct[match("LesionC", lev_les)], 4),
    pct_min     = round(min(pct), 4),
    pct_max     = round(max(pct), 4),
    max_minus_min_pct = round(max(pct) - min(pct), 4),
    max_fold    = round(max(pct) / max(min(pct), 1e-12), 3),
    test        = test,
    p_value     = p,
    cramers_V   = round(V, 4),
    min_expected = round(exp_min, 3),
    stringsAsFactors = FALSE
  )
}
st <- do.call(rbind, res)
st$FDR <- round(p.adjust(st$p_value, method = "BH"), 6)
st <- st[order(-st$n_total), ]
write.csv(st, file.path(out_fig, "SC_celltype_by_lesion_stats.csv"), row.names = FALSE)
cat("\n[written] SC_celltype_by_lesion_stats.csv :", nrow(st), "x", ncol(st), "\n")
print(st, row.names = FALSE)

# ---- T3: locate Treg / plasma / mast in celltype vocabulary ----
cat("\n=== T3: celltype vocabulary search ===\n")
vocab <- sort(unique(md_c$celltype))
cat("all celltypes:\n"); print(vocab)
pat <- list(
  Treg   = c("treg", "regulatory", "T_reg", "Tregs"),
  plasma = c("plasma", "plasmablast", "Plasma_cell"),
  mast   = c("mast")
)
for (nm in names(pat)) {
  hit <- vocab[Reduce(`|`, lapply(pat[[nm]], function(p) grepl(p, vocab, ignore.case = TRUE)))]
  cat(sprintf("%-7s exact-vocabulary match: %s\n", nm,
              if (length(hit) == 0) "NONE" else paste(hit, collapse = ", ")))
}

# closest available coarse proxy
proxy <- list(
  Treg   = grep("^T_cell", vocab, value = TRUE, ignore.case = TRUE),
  plasma = grep("B_cell|Plasma", vocab, value = TRUE, ignore.case = TRUE),
  mast   = "NONE"
)
for (nm in names(proxy)) {
  cat(sprintf("%-7s coarse proxy: %s\n", nm,
              if (length(proxy[[nm]]) == 0) "NONE" else paste(proxy[[nm]], collapse = ", ")))
}

# cluster-level marker check for the three lineages
mk_file <- file.path(RES, "sc_cluster_markers.csv")
cat("\n=== cluster-level marker evidence (from sc_cluster_markers.csv) ===\n")
if (file.exists(mk_file)) {
  mk <- read.csv(mk_file, stringsAsFactors = FALSE)
  keyg <- c("FOXP3", "IL2RA", "CTLA4", "IKZF2",   # Treg
            "MZB1", "JCHAIN", "SDC1", "XBP1", "IGHG1",  # plasma
            "TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2") # mast
  sub <- mk[mk$gene %in% keyg, c("cluster", "gene", "avg_log2FC", "p_val_adj")]
  if (nrow(sub) > 0) {
    sub <- sub[order(sub$cluster, sub$gene), ]
    print(sub, row.names = FALSE)
  } else cat("none of the key markers appear in the stored marker table\n")
} else cat("marker file missing\n")

# ---- T4: keratinocyte / T-cell share per lesion ----
cat("\n=== T4: keratinocyte-dominance check ===\n")
ker <- intersect(c("Keratinocytes", "Epithelial_cells"), vocab)
tcl <- intersect(c("T_cells", "T_cell", "NK_cell"), vocab)
cat("keratinocyte-lineage labels used:", paste(ker, collapse = ", "), "\n")
cat("T-cell labels used:", paste(intersect("T_cells", vocab), collapse = ", "), "; T+NK =", paste(tcl, collapse = ", "), "\n")
for (nm in c("keratinocyte_lineage", "T_cells_only", "T_plus_NK")) {
  ids <- switch(nm, keratinocyte_lineage = ker, T_cells_only = intersect("T_cells", vocab), T_plus_NK = tcl)
  if (length(ids) == 0) next
  v <- colSums(ct[ids, , drop = FALSE])
  cat(sprintf("%-22s n: A=%d B=%d C=%d | pct: A=%.2f B=%.2f C=%.2f | overall=%.2f%%\n",
              nm, v[1], v[2], v[3],
              100*v[1]/N_les[1], 100*v[2]/N_les[2], 100*v[3]/N_les[3],
              100*sum(v)/sum(N_les)))
}

cat("\n=== DONE ===\n")
sessionInfo()
