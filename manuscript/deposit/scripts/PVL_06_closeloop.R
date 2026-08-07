suppressMessages({
  library(Seurat); library(data.table); library(ggplot2); library(pheatmap)
  library(RColorBrewer)
})
## 输入数据目录 = .env 的 DATASET_PATH（WSL 下为 /mnt/g/...）；分析结果仍留本地 ~/PVL/results
DATASET_PATH <- "/mnt/g/CMCC Sync/Medical Data/Stomatology/PVL"
RAW   <- file.path(DATASET_PATH, "raw")
RES   <- "~/PVL/results"
MET   <- file.path(DATASET_PATH, "metadata")

## ---------- load discovery signature (GSE26549, 101 genes) ----------
deg  <- fread(file.path(RES, "GSE26549_progression_DEG.csv"))
sig  <- unique(deg[adj.P.Val < 0.05 & abs(logFC) > 0.5]$gene); sig <- sig[!is.na(sig)]
w    <- setNames(deg$logFC, deg$gene)[sig]; w <- w[!is.na(w)]
cat("discovery signature genes:", length(sig), "\n")

## ---------- load single-cell object ----------
rds  <- file.path(RES, "PVL_sc_combined.rds")
so   <- readRDS(rds)
cat("sc object:", ncol(so), "cells x", nrow(so), "genes; layers:",
    paste(SeuratObject::Layers(so), collapse = "/"), "\n")
if (!"data" %in% SeuratObject::Layers(so)) so <- JoinLayers(so)
feat <- rownames(so)
sig_in <- intersect(sig, feat)
cat("signature genes present in sc object:", length(sig_in), "/", length(sig), "\n")

## (1) per-cell direct z-score signature score (single-cell analog of bulk z-score)
mat    <- as.matrix(LayerData(so, layer = "data")[sig_in, ])   # genes x cells, log-norm
sig_z  <- t(scale(t(mat)))                                     # z per gene
cell_score <- colMeans(sig_z, na.rm = TRUE)
so$ProgSig <- cell_score
df <- data.frame(cell = colnames(so), lesion = so$lesion, score = so$ProgSig,
                 cluster = as.character(Idents(so)))
lvl <- c("LesionA","LesionB","LesionC")
df$lesion <- factor(df$lesion, levels = lvl)
m <- aggregate(score ~ lesion, df, mean)$score
names(m) <- lvl
cat("per-cell signature z-score mean by lesion:\n"); print(round(m,4))
# effect size: eta^2 of lesion on score + KW
kw  <- kruskal.test(score ~ lesion, df)
g   <- split(df$score, df$lesion)
grand <- mean(df$score); ssb <- sum(sapply(g, function(x) length(x)*(mean(x)-grand)^2))
sst <- sum((df$score-grand)^2); eta2 <- ssb/sst
cat("KW p =", format.pval(kw$p.value, digits=3), "  lesion eta^2 (effect size) =", round(eta2,4), "\n")
p1 <- ggplot(df, aes(x = lesion, y = score, fill = lesion)) +
  geom_violin(scale = "width", alpha = .6) + geom_boxplot(width = .15, outlier.size = .3) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  labs(title = "Bulk progression-signature score per cell (z per gene)",
       subtitle = paste0("lesion eta^2 = ", round(eta2,3), " (negligible); not a within-lesion gradient"),
       y = "mean z-score across 87 signature genes", x = "") +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(RES, "PVL6_sigscore_by_lesion.png"), p1, width = 5, height = 5)

## identify signature-high clusters
Idents(so) <- "seurat_clusters"
cl_score <- aggregate(score ~ cluster, df, mean)
cl_score$seurat_clusters <- cl_score$cluster
top_cl <- cl_score[order(-cl_score$score), ][1:5, ]
cat("top 5 signature-high clusters:\n"); print(data.frame(cluster = top_cl$seurat_clusters, score = round(top_cl$score, 4)))

## (2) signature-gene expression heatmap across 24 clusters ----------
av <- as.matrix(AverageExpression(so, features = sig_in, group.by = "seurat_clusters")$RNA)
# strip Seurat 'g' prefix on numeric cluster ids
cc <- sub("^g", "", colnames(av)); colnames(av) <- cc
cc <- cc[order(as.integer(cc))]; av <- av[, cc]
av_scaled <- t(scale(t(av)))   # z per gene
# dominant lesion per cluster (correct table column name)
comp <- as.data.frame(table(so$seurat_clusters, so$lesion))
# dominant lesion per cluster
dom <- lapply(split(comp, comp$Var1), function(d) as.character(d$Var2[which.max(d$Freq)]))
cl2les <- unlist(dom)
domles <- cl2les[colnames(av_scaled)]
ann_col <- data.frame(Lesion = factor(domles, levels = lvl))
rownames(ann_col) <- colnames(av_scaled)
pal <- brewer.pal(3, "Set2"); names(pal) <- lvl
pheatmap(av_scaled, cluster_cols = FALSE, cluster_rows = TRUE,
         show_rownames = nrow(av_scaled) <= 120,
         annotation_col = ann_col, annotation_colors = list(Lesion = pal),
         fontsize = 7, main = "Signature-gene mean expression per cluster (z per gene)",
         filename = file.path(RES, "PVL6_sig_heatmap_clusters.png"), width = 7, height = 9)
cat("heatmap: clusters =", ncol(av_scaled), " signature genes plotted =", nrow(av_scaled), "\n")

## (3) batch-effect evidence: lesion PCA variance + split UMAP ----------
set.seed(1)
so <- RunPCA(so, features = VariableFeatures(so), npcs = 20, verbose = FALSE)
emb <- Embeddings(so, "pca")[, 1:10]
eta2pc <- sapply(1:ncol(emb), function(k){
  g <- split(emb[, k], so$lesion)
  grand <- mean(emb[, k]); ssb <- sum(sapply(g, function(x) length(x)*(mean(x)-grand)^2))
  sst <- sum((emb[, k]-grand)^2); ssb/sst
})
cat("lesion eta^2 on PC1..PC10:\n"); print(round(eta2pc,3))
cat("mean lesion eta^2 (top10 PCs):", round(mean(eta2pc),3), "\n")
write.csv(data.frame(PC = paste0("PC",1:10), lesion_eta2 = eta2pc),
          file.path(RES, "PVL6_lesion_pca_variance.csv"), row.names = FALSE)
p3 <- DimPlot(so, reduction = "umap", group.by = "lesion", split.by = "lesion",
              cols = pal, ncol = 3) + theme(legend.position = "none") +
  labs(title = "UMAP split by lesion (batch ~ biology)")
ggsave(file.path(RES, "PVL6_split_umap_lesion.png"), p3, width = 11, height = 4)

## (4) LASSO-20 vs cluster-marker overlap ----------
lasso <- fread(file.path(RES, "PVL4_lasso_coefficients.csv"))
lasso_genes <- lasso$gene[!is.na(lasso$gene) & lasso$gene != ""]
markers <- fread(file.path(RES, "sc_cluster_markers.csv"))
ov <- intersect(lasso_genes, unique(markers$gene))
cat("LASSO genes:", length(lasso_genes), "; overlap with sc markers:", length(ov),
    " ->", if(length(ov)) paste(ov, collapse = ",") else "NONE", "\n")
if (length(ov)) {
  sub <- markers[gene %in% ov, .(gene, cluster, avg_log2FC)]
  best <- sub[, .SD[which.max(avg_log2FC)], by = gene]
  write.csv(best, file.path(RES, "PVL6_lasso_marker_overlap.csv"), row.names = FALSE)
} else {
  write.csv(data.frame(gene=character(), cluster=integer(), avg_log2FC=numeric()),
            file.path(RES, "PVL6_lasso_marker_overlap.csv"), row.names = FALSE)
}

## ---------- summary ----------
summ <- data.frame(
  metric = c("sig_genes_in_sc","score_LesionA","score_LesionB","score_LesionC",
             "KW_p","lesion_eta2","top_cluster_score","top_cluster_id",
             "lesion_eta2_mean_PC1_10","LASSO_marker_overlap_n"),
  value = c(length(sig_in), round(m["LesionA"],4), round(m["LesionB"],4), round(m["LesionC"],4),
            round(kw$p.value,4), round(eta2,4),
            round(top_cl$score[1],4), as.character(top_cl$seurat_clusters[1]),
            round(mean(eta2pc),3), length(ov)))
write.csv(summ, file.path(RES, "PVL6_summary.csv"), row.names = FALSE)
cat("DONE Phase6\n")
