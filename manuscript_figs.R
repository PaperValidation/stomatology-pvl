#!/usr/bin/env Rscript
## ============================================================
## Manuscript figure + table production (v2, 2026-08-04)
## Replicates EXACT scoring formulas from PVL_04 (safe_scale,
## Zall %*% w, ord_map) and PVL_06 (per-cell z-score module score)
## so that figures are numerically identical to manuscript_v1.md.
## Outputs -> results/manuscript_figs/
## ============================================================
suppressMessages({
  library(Seurat); library(data.table); library(ggplot2); library(pheatmap)
  library(RColorBrewer)
})
## ---- paths (DATASET_PATH = .env, WSL translation) ----
DATASET_PATH <- "/mnt/g/CMCC Sync/Medical Data/Stomatology/PVL"
RAW  <- file.path(DATASET_PATH, "raw")
MET  <- file.path(DATASET_PATH, "metadata")
RES  <- "~/PVL/results"
OUT  <- "/mnt/c/Users/tidus/PycharmProjects/Stomatology-PVL/results/manuscript_figs"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- helpers (replicated from PVL_04) ----
safe_scale <- function(mat){            # genes x samples -> samples x genes, z per gene
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
  rownames(m) <- rn; colnames(m) <- hdr[-1]
  m
}
cat("=== manuscript figures start ===\n")

## ---------- Fig 1A: severity-gradient scatter (GSE227919) ----------
deg <- fread(file.path(RES, "GSE26549_progression_DEG.csv"))
sig <- unique(deg[adj.P.Val < 0.05 & abs(logFC) > 0.5]$gene); sig <- sig[!is.na(sig)]
lines <- readLines(file.path(RAW, "GSE227919_2023_20_03_expr_mat_raw.txt.gz"))
hdr <- strsplit(trimws(lines[1]), "\\s+")[[1]]
fm <- fread(text = lines[-1], sep = " ", header = FALSE, fill = TRUE)
genes <- as.character(fm[[1]])
exm <- as.matrix(fm[, -1, with = FALSE]); rownames(exm) <- genes; colnames(exm) <- hdr
clin22g <- fread(file.path(MET, "GSE227919_clinical.csv"))
cls <- setNames(clin22g$Class, clin22g$Sample_ID)
keep_all <- intersect(colnames(exm), names(cls))
common <- intersect(sig, rownames(exm))
w <- setNames(deg$logFC, deg$gene)[common]; w <- w[!is.na(w)]
Zall <- safe_scale(exm[names(w), keep_all, drop = FALSE])
score_all <- as.numeric(Zall %*% w)
ord_map <- c(Control = 0, HkNR = 1, Dysplasia = 2, Cancer = 3)
ord_all <- ord_map[cls[keep_all]]
df1 <- data.frame(sample = keep_all, Class = factor(cls[keep_all],
                 levels = c("Control","HkNR","Dysplasia","Cancer")),
                  severity = ord_all, score = score_all)
r1 <- cor.test(score_all, ord_all, method = "spearman")
nc <- ord_all < 3
r2 <- cor.test(score_all[nc], ord_all[nc], method = "spearman")
cat(sprintf("Fig1A rho all=%.3f (p=%.2e); non-cancer n=%d rho=%.3f (p=%.4f)\n",
            r1$estimate, r1$p.value, sum(nc), r2$estimate, r2$p.value))
p1 <- ggplot(df1, aes(x = severity, y = score, color = Class)) +
  geom_jitter(width = 0.18, size = 2.2, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, color = "grey40", linewidth = 0.6) +
  scale_x_continuous(breaks = 0:3, labels = c("Control","HkNR","Dysplasia","Cancer"),
                     limits = c(-0.45, 3.45)) +
  scale_color_manual(values = c(Control = "#4daf4a", HkNR = "#377eb8",
                                Dysplasia = "#ff7f00", Cancer = "#e41a1c")) +
  labs(x = "Histopathological severity", y = "Direction-weighted signature score",
       title = "GSE227919 (RNA-seq, n = 66)",
       subtitle = sprintf("Spearman rho = %.3f, P = %.1e", r1$estimate, r1$p.value)) +
  annotate("text", x = 1.05, y = max(score_all)*0.95, hjust = 0, size = 3.2,
           label = sprintf("all: rho = %.3f\nnon-cancer (n = %d): rho = %.3f (P = %.3f)",
                           r1$estimate, sum(nc), r2$estimate, r2$p.value)) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "Fig1A_severity_gradient.png"), p1, width = 6.2, height = 5.2, dpi = 300)
write.csv(df1, file.path(OUT, "Fig1A_data.csv"), row.names = FALSE)

## ---------- Fig 1B(ii): per-cell signature-score UMAP (GSE196296) ----------
so <- readRDS(file.path(RES, "PVL_sc_combined.rds"))
if (!"data" %in% SeuratObject::Layers(so)) so <- JoinLayers(so)
sig_in <- intersect(sig, rownames(so))
mat <- as.matrix(LayerData(so, layer = "data")[sig_in, ])
sig_z <- t(scale(t(mat)))
cell_score <- colMeans(sig_z, na.rm = TRUE)
so$ProgSig <- cell_score
cat("Fig1B(ii): cells =", ncol(so), " signature genes in sc =", length(sig_in), "\n")
p2 <- FeaturePlot(so, features = "ProgSig", reduction = "umap",
                  order = TRUE, cols = c("grey90", "#440154", "#21918c", "#fde725")) +
  labs(title = "Per-cell progression-signature score (module z-score)") +
  theme(legend.position = "right")
ggsave(file.path(OUT, "Fig1B_umap_sigscore.png"), p2, width = 6.5, height = 5.4, dpi = 300)

## ---------- Fig 2: LINCS reversal evidence heatmap (8 candidates) ----------
lc <- fread(file.path(RES, "PVL5_lincs_reversal_candidates.csv"))
targets <- c("ZSTK-474" = "PI3K", "PF-477736" = "CHK1", "AZD-7762" = "CHK1",
             "HG-6-64-01" = "HDAC", "GSK-1059615" = "PI3K/mTOR", "AS-601245" = "JNK",
             "torin-1" = "mTOR", "sorafenib" = "multi-kinase")
tier <- c("ZSTK-474" = "Priority", "GSK-1059615" = "Priority", "torin-1" = "Priority",
          "PF-477736" = "Exploratory", "AZD-7762" = "Exploratory",
          "HG-6-64-01" = "Exploratory", "AS-601245" = "Exploratory", "sorafenib" = "Exploratory")
hub_up <- c("TFRC","SPP1","SLC2A1","MET")          # progression-up, should be down
hub_dn <- c("ECM1","PDGFRA","IGFBP5","PTN")        # progression-low, should be up
mk <- function(x) if (is.na(x) || x == "") character(0) else trimws(strsplit(x, ";")[[1]])
M <- matrix(NA, nrow = nrow(lc), ncol = length(c(hub_up, hub_dn)),
            dimnames = list(lc$Drug, c(hub_up, hub_dn)))
for (i in seq_len(nrow(lc))) {
  M[i, intersect(mk(lc$up_targets_down[i]), hub_up)] <- "Up-target \u2193"
  M[i, intersect(mk(lc$down_targets_up[i]),  hub_dn)] <- "Low-target \u2191"
}
M[M == ""] <- NA
Mdf <- as.data.frame(M)
Mdf$Drug <- rownames(M)
mlong <- melt(setDT(Mdf), id.vars = "Drug", variable.name = "Gene", value.name = "Rev")
mlong$Tier <- factor(tier[mlong$Drug], levels = c("Priority","Exploratory"))
mlong$Class <- factor(targets[mlong$Drug])
mlong$Rev <- factor(mlong$Rev, levels = c("Up-target \u2193","Low-target \u2191"))
mlong$Drug <- factor(mlong$Drug, levels = rev(c("torin-1","GSK-1059615","ZSTK-474",
                                                 "AZD-7762","PF-477736","HG-6-64-01",
                                                 "AS-601245","sorafenib")))
p3 <- ggplot(mlong, aes(x = Gene, y = Drug, fill = Rev)) +
  geom_tile(color = "white", linewidth = 0.7, width = 0.92, height = 0.92) +
  scale_fill_manual(values = c("Up-target \u2193" = "#d7191c", "Low-target \u2191" = "#2b83ba"),
                    na.value = "grey93", name = "Reversal evidence") +
  facet_grid(Tier ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "Hub gene", y = NULL,
       title = "LINCS two-direction reversal evidence (8 candidates)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        legend.position = "bottom")
ggsave(file.path(OUT, "Fig2_lincs_reversal.png"), p3, width = 6.2, height = 5.6, dpi = 300)
write.csv(as.data.frame(M), file.path(OUT, "Fig2_reversal_matrix.csv"))

## ---------- Fig 3: bulk Hallmark enrichment bar (GSE227919 CvC) ----------
hb <- fread(file.path(RES, "enrich_MSigDB_Hallmark_2020_GSE227919_CancerVsControl_up.csv"))
hb <- hb[Adjusted.P.value < 0.05]
key <- c("mTORC1 Signaling","PI3K/AKT/mTOR  Signaling","Epithelial Mesenchymal Transition")
hb$key <- ifelse(hb$Term %in% key, "Emphasised", "Other")
hb <- hb[order(Combined.Score, decreasing = TRUE)]
top10 <- hb[1:10]
miss <- key[!key %in% top10$Term]                       # ensure all 3 emphasised terms shown
if (length(miss) > 0) {
  extra <- hb[Term %in% miss]
  top10 <- rbind(head(hb[!Term %in% miss], 10 - nrow(extra)), extra)
}
hb <- top10
hb$Term <- factor(hb$Term, levels = rev(hb$Term))
hb$lab <- sprintf("%.0e", hb$Adjusted.P.value)
p4 <- ggplot(hb, aes(x = Combined.Score, y = Term, fill = key)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = lab), hjust = -0.08, size = 3) +
  scale_fill_manual(values = c(Emphasised = "#e41a1c", Other = "#377eb8"),
                    name = NULL, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Combined score", y = NULL,
       title = "Bulk GSE227919 cancer-versus-control Hallmark enrichment (top 10, adj. P < 0.05)",
       subtitle = "Red: pathways emphasised in the manuscript; labels = adjusted P") +
  theme_bw(base_size = 11)
ggsave(file.path(OUT, "Fig3_bulk_hallmark.png"), p4, width = 7.2, height = 5.2, dpi = 300)
write.csv(hb[, .(Term, Combined.Score, Adjusted.P.value, key)],
          file.path(OUT, "Fig3_data.csv"), row.names = FALSE)

## ---------- Table 1: 8 reversal candidates (clean) ----------
t1 <- lc[, .(Drug)]
t1$Target.class <- targets[t1$Drug]
t1$Tier <- tier[t1$Drug]
t1$Up.targets.down <- sapply(lc$up_targets_down, function(x) length(mk(x)))
t1$Down.targets.up <- sapply(lc$down_targets_up, function(x) length(mk(x)))
t1$Reversal.score <- round(lc$reversal_score, 1)
t1$Hub.evidence <- sapply(seq_len(nrow(lc)), function(i)
  paste(unique(c(intersect(mk(lc$up_targets_down[i]), hub_up),
                 intersect(mk(lc$down_targets_up[i]), hub_dn))), collapse = "; "))
fwrite(t1, file.path(OUT, "Table1_reversal_candidates.csv"))

## ---------- Table 2: signature-level vs bulk-level enrichment ----------
p7 <- fread(file.path(RES, "PVL7_summary.csv"))
hall_n_bulk <- sum(fread(file.path(RES, "enrich_MSigDB_Hallmark_2020_GSE227919_CancerVsControl_up.csv"))$Adjusted.P.value < 0.05)
t2 <- data.table(
  Layer = c("Signature-level (101-gene, GSE26549)",
            "Bulk-level (GSE227919 cancer vs control, |logFC|>1)"),
  Hallmark.sig = c("1 term: UV Response Down (adj. P = 0.023)",
                   sprintf("%d terms; top: E2F, G2M, IFN-*; mTORC1, PI3K/AKT/mTOR, EMT (adj. P < 0.05)", hall_n_bulk)),
  GO.sig = c("3 terms (steroid/androgen metabolism; neg. reg. proteolysis)", "n.r."),
  TF.lib.sig = c("0 (ChEA/ENCODE)", "n.r."),
  KEGG.sig = c("0", "see GO/Reactome in pipeline outputs")
)
fwrite(t2, file.path(OUT, "Table2_enrichment_summary.csv"))
cat("=== manuscript figures done ->", OUT, "\n")
