#!/usr/bin/env Rscript
## ============================================================
## PVL Phase 3: single-cell RNA-seq (GSE196296, 10x)
## 3 samples = same PVL lesion pathological continuum
##   LesionA (KUS) -> LesionB (corrugated) -> LesionC (bulky hyperkeratosis)
## ============================================================
Sys.setenv(http_proxy  = "socks5://127.0.0.1:8428",
           https_proxy = "socks5://127.0.0.1:8428",
           ftp_proxy   = "socks5://127.0.0.1:8428")
options(download.file.method = "libcurl")

suppressMessages({
  library(Seurat); library(Matrix); library(data.table); library(dplyr)
  library(ggplot2); library(cowplot)
})

## 输入数据目录 = .env 的 DATASET_PATH（WSL 下为 /mnt/g/...）；分析结果仍留本地 ~/PVL/results
DATASET_PATH <- "/mnt/g/CMCC Sync/Medical Data/Stomatology/PVL"
RAW  <- file.path(DATASET_PATH, "raw")
RES  <- "~/PVL/results"
dir.create(RES, showWarnings = FALSE)
cat("=== PVL Phase 3 start ===\n")

samples <- list(
  LesionA = "GSM5870060",
  LesionB = "GSM5870062",
  LesionC = "GSM5870063"
)

read_10x <- function(gsm, lesion){
  matp  <- file.path(RAW, sprintf("%s_%s_matrix.mtx.gz",  gsm, lesion))
  bcp   <- file.path(RAW, sprintf("%s_%s_barcodes.tsv.gz", gsm, lesion))
  ftp   <- file.path(RAW, sprintf("%s_%s_features.tsv.gz", gsm, lesion))
  ft    <- read.delim(ftp, header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  sym   <- ifelse(is.na(ft$V2) | ft$V2 == "", ft$V1, ft$V2)
  sym   <- make.unique(sym)
  mat   <- readMM(gzfile(matp))                      # features x cells
  rownames(mat) <- sym
  bc    <- readLines(gzfile(bcp))
  colnames(mat) <- paste0(bc, "_", lesion)
  so    <- CreateSeuratObject(counts = mat, project = lesion)
  so$lesion <- rep(lesion, ncol(so))
  cat(sprintf("  %s: %d genes x %d cells\n", lesion, nrow(so), ncol(so)))
  so
}

objs <- lapply(names(samples), function(les) read_10x(samples[[les]], les))
names(objs) <- names(samples)
combined <- merge(objs[[1]], objs[2:3], add.cell.ids = names(samples))
cat("merged:", nrow(combined), "genes x", ncol(combined), "cells\n")

## ---------- QC ----------
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 25)
cat("after QC:", ncol(combined), "cells\n")

## ---------- standard workflow ----------
combined <- NormalizeData(combined, verbose = FALSE) %>%
  FindVariableFeatures(nfeatures = 3000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 50, verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.6, verbose = FALSE)
combined <- JoinLayers(combined)   # consolidate per-sample layers (data.LesionA/B/C -> single data)
cat("clusters:", length(unique(combined$seurat_clusters)), "\n")

## ---------- lesion composition per cluster ----------
comp <- table(combined$seurat_clusters, combined$lesion)
write.csv(as.data.frame(comp), file.path(RES, "sc_cluster_lesion_composition.csv"))
cat("lesion composition per cluster:\n"); print(comp)

## ---------- SingleR annotation (if celldex available) ----------
combined$celltype <- "unassigned"
if (suppressWarnings(requireNamespace("SingleR", quietly = TRUE)) &&
    suppressWarnings(requireNamespace("celldex", quietly = TRUE))) {
  cat("running SingleR annotation...\n")
  ref <- tryCatch(celldex::HumanPrimaryCellAtlasData(), error = function(e){ cat("ref dl err:", e$message, "\n"); NULL })
  if (!is.null(ref)) {
    test.mat <- LayerData(combined, layer = "data")   # sparse; SingleR handles dgCMatrix
    pred <- SingleR::SingleR(test = test.mat, ref = ref, labels = ref$label.main)
    combined$celltype <- pred$pruned.labels
    write.csv(as.data.frame(table(combined$celltype)),
              file.path(RES, "sc_celltype_counts.csv"), row.names = FALSE)
    cat("cell types assigned.\n")
  }
} else {
  cat("celldex/SingleR unavailable -> skipping annotation\n")
}

## early save: persist the expensive object before the slow marker step
saveRDS(combined, file.path(RES, "PVL_sc_combined.rds"))
cat("object saved (early).\n")

## ---------- outputs ----------
p1 <- DimPlot(combined, group.by = "lesion",   pt.size = 0.3) + ggtitle("UMAP by lesion (pathology continuum)")
p2 <- DimPlot(combined, group.by = "seurat_clusters", pt.size = 0.3) + ggtitle("UMAP by cluster")
plots <- list(p1, p2)
if (any(combined$celltype != "unassigned")) {
  p3 <- DimPlot(combined, group.by = "celltype", pt.size = 0.3, label = TRUE) + ggtitle("UMAP by cell type (SingleR)")
  plots <- c(plots, list(p3))
}
pdf(file.path(RES, "sc_umap.pdf"), width = 10, height = 5)
  print(plot_grid(plotlist = plots, nrow = 1))
dev.off()
png(file.path(RES, "sc_umap.png"), width = 1600, height = 800, res = 120)
  print(plot_grid(plotlist = plots, nrow = 1))
dev.off()

## cluster markers (restricted to variable features + downsampled for speed)
markers <- FindAllMarkers(combined, features = VariableFeatures(combined),
                          max.cells.per.ident = 200, only.pos = TRUE,
                          min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
write.csv(markers, file.path(RES, "sc_cluster_markers.csv"), row.names = FALSE)
cat("markers written:", nrow(markers), "\n")

## save object
saveRDS(combined, file.path(RES, "PVL_sc_combined.rds"))
cat("\n=== PVL Phase 3 DONE ===\n")
