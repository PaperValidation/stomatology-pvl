#!/usr/bin/env Rscript
## ============================================================
## PVL Phase 10: quantification of single-cell signature localisation
## For each seurat cluster, test the direction-weighted signature score
## (cluster vs rest) with Wilcoxon rank-sum + BH-FDR, and report an
## ordinal effect size (Cliff's delta) and binary-membership eta^2.
## Reads only local files (no network).
## ============================================================
suppressMessages({
  library(Seurat); library(data.table); library(ggplot2)
})

DATASET_PATH <- "/mnt/g/CMCC Sync/Medical Data/Stomatology/PVL"
RAW   <- file.path(DATASET_PATH, "raw")
RES   <- "~/PVL/results"

## ---------- load discovery signature (GSE26549, 101 genes) ----------
deg  <- fread(file.path(RES, "GSE26549_progression_DEG.csv"))
sig  <- unique(deg[adj.P.Val < 0.05 & abs(logFC) > 0.5]$gene); sig <- sig[!is.na(sig)]
w    <- setNames(deg$logFC, deg$gene)[sig]; w <- w[!is.na(w)]

## ---------- load single-cell object ----------
so <- readRDS(file.path(RES, "PVL_sc_combined.rds"))
if (!"data" %in% SeuratObject::Layers(so)) so <- JoinLayers(so)
feat   <- rownames(so)
sig_in <- intersect(sig, feat)
cat("signature genes present in sc object:", length(sig_in), "/", length(sig), "\n")

## ---------- per-cell direction-weighted z-score (same as PVL_06) ----------
mat     <- as.matrix(LayerData(so, layer = "data")[sig_in, ])
sig_z   <- t(scale(t(mat)))
cell_score <- colMeans(sig_z, na.rm = TRUE)
so$ProgSig <- cell_score

df <- data.frame(cell   = colnames(so),
                 lesion = as.character(so$lesion),
                 cluster = as.character(Idents(so)),
                 score  = so$ProgSig, stringsAsFactors = FALSE)

## ---------- per-cluster statistics: cluster vs rest ----------
clusters_int <- sort(unique(df$cluster), decreasing = FALSE)
rows <- list()
for (cl in clusters_int) {
  ingrp <- df$score[df$cluster == cl]
  outgrp <- df$score[df$cluster != cl]
  nx <- length(ingrp); ny <- length(outgrp)
  wt <- wilcox.test(ingrp, outgrp, exact = FALSE)
  ## Cliff's delta = P(ingrp>outgrp) - P(outgrp>ingrp)  (ties ignored, small bias)
  R_in <- sum(rank(c(ingrp, outgrp))[1:nx])
  U_in <- R_in - nx * (nx + 1) / 2
  P_gt <- U_in / (nx * ny)
  cliffs <- 2 * P_gt - 1
  ## eta^2 of binary membership on score
  grand <- mean(df$score)
  ssb <- nx * (mean(ingrp) - grand)^2 + ny * (mean(outgrp) - grand)^2
  sst <- sum((df$score - grand)^2)
  eta2 <- ssb / sst
  rows[[cl]] <- data.frame(cluster = cl, n_cells = nx,
                           mean_score = round(mean(ingrp), 4),
                           wilcox_p = wt$p.value,
                           cliffs_delta = round(cliffs, 4),
                           eta2 = round(eta2, 5),
                           stringsAsFactors = FALSE)
}
stat_df <- rbindlist(rows)
stat_df$FDR <- p.adjust(stat_df$wilcox_p, method = "BH")
stat_df <- stat_df[order(-stat_df$mean_score)]
fwrite(stat_df, file.path(RES, "PVL10_cluster_sigscore_stats.csv"), row.names = FALSE)

cat("=== top 6 signature-high clusters (cluster vs rest) ===\n")
print(stat_df[1:6, .(cluster, n_cells, mean_score, wilcox_p, FDR, cliffs_delta, eta2)])

nsig <- sum(stat_df$FDR < 0.05)
cat("\nclusters significant at FDR<0.05:", nsig, "of", nrow(stat_df), "\n")

named <- c("5", "13", "22", "23", "15")
cat("\n=== the 5 previously named high clusters ===\n")
print(stat_df[cluster %in% named, .(cluster, n_cells, mean_score, wilcox_p, FDR, cliffs_delta, eta2)])

## ---------- bar plot: mean score per cluster, significance-coloured ----------
p <- ggplot(stat_df, aes(x = reorder(cluster, mean_score), y = mean_score, fill = FDR < 0.05)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("grey70", "steelblue"),
                    labels = c("FDR>=0.05", "FDR<0.05")) +
  labs(title = "Per-cluster progression-signature score (direction-weighted z-score)",
       subtitle = "blue = significant vs remaining cells (Wilcoxon BH-FDR<0.05)",
       x = "seurat cluster", y = "mean per-cell signature score") +
  theme_minimal() + theme(legend.title = element_blank())
ggsave(file.path(RES, "PVL10_cluster_sigscore_bar.png"), p, width = 6, height = 7)

## ---------- summary for manuscript ----------
nsig_clusters <- stat_df[FDR < 0.05, cluster]
summ <- data.frame(
  metric = c("n_clusters", "n_cells_total", "n_clusters_FDR0.05",
             "top_cluster", "top_cluster_mean", "top_cluster_cliffs",
             "top_cluster_FDR",
             "named5_FDR0.05_count", "named5_clusters"),
  value = c(length(clusters_int), nrow(df), nsig,
            stat_df$cluster[1], stat_df$mean_score[1], stat_df$cliffs_delta[1],
            format.pval(stat_df$FDR[1], digits = 3),
            sum(stat_df[cluster %in% named]$FDR < 0.05),
            paste(stat_df[cluster %in% named]$cluster, collapse = ",")))
fwrite(summ, file.path(RES, "PVL10_summary.csv"), row.names = FALSE)
cat("\nDONE PVL10\n")
